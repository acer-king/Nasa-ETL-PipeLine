#!/usr/bin/env python3
import unittest
from random import randint
import pyutils
from pyutils import re_snowflake
import pandas as pd


class TestSnowflake(unittest.TestCase):
    def setUp(self) -> None:
        CREDENTIALS = pyutils.get_secret('redas/prod')
        self.sf = re_snowflake.ReSnowflake(warehouse="DATALOADING",
                                           database="TESTDB",
                                           schema="PUBLIC",
                                           role='READ_WRITE_RL',
                                           sf_user=CREDENTIALS['SNOWFLAKE_USER'],
                                           sf_password=CREDENTIALS['SNOWFLAKE_PASS'],
                                           sf_account=CREDENTIALS['SNOWFLAKE_ACCOUNT'])

    def test_stage_upsert(self):
        test_table = "TESTTABLE{}".format(str(randint(100, 999)))
        create_test_table_query = 'create or replace table {} ("A" TEXT, "B" NUMERIC, ' \
                                  'PRIMARY KEY("A"));'.format(test_table)
        cs = self.sf.con.cursor()
        try:
            cs.execute(create_test_table_query)
        except Exception as e:
            self.sf.logging.exception("Error in creating test table: \n {}".format(e))
            return False
        finally:
            cs.close()
        # Test inserting data into an empty table using stage_upsert.
        df1 = pd.DataFrame({'A': ['a', 'b'], 'B': [1, 2]})
        self.sf.stage_upsert(data=df1, table=test_table)

        retrieve_df1 = self.sf.read_sql('SELECT * FROM {};'.format(test_table))
        pd.testing.assert_frame_equal(retrieve_df1, df1)

        # Test if data is updated using upsert.
        df2 = pd.DataFrame({'A': ['b', 'c'], 'B': [3, 4]})
        self.sf.stage_upsert(data=df2, table=test_table)

        expected_df = pd.DataFrame({'A': ['a', 'b', 'c'], 'B': [1, 3, 4]})
        retrieve_df2 = self.sf.read_sql('SELECT * FROM {} ORDER BY A;'.format(test_table))
        pd.testing.assert_frame_equal(expected_df, retrieve_df2)

