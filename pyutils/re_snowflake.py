from . import *
import pandas as pd
import snowflake.connector
from snowflake.connector import DictCursor
from typing import List
from random import randint
from pathlib import Path
import tempfile

# initialize logging
logging.basicConfig(format='%(asctime)s  %(levelname)s: %(message)s', level=logging.INFO)


class ReSnowflake:

    def __init__(self, warehouse, database, role, schema, sf_user, sf_password, sf_account):
        # create connection start warehouse and connect to DB/schema
        self.con = None
        self.con = snowflake.connector.connect(
            user=sf_user,
            password=sf_password,
            account=sf_account,
            region='us-east-1',
            warehouse=warehouse,
            database=database,
            role=role,
            schema=schema,
            timezone='UTC'
        )

    def __del__(self):
        if self.con:
            del self.con

    def read_sql(self, query: str) -> pd.DataFrame:
        """
        Return results of a SQL query to snowflake as a DataFrame.
        -----------------------------------------------------------
        :param query: String
                        SQL query
        :return: pandas.DataFrame
        """

        cs = self.con.cursor(DictCursor)
        try:
            cs.execute(query)
            data = cs.fetchall()
            df = pd.DataFrame(data)
        except Exception as e:
            logging.exception("Error in loading query: \n {}".format(e))
            cs.close()
            return 1
        finally:
            cs.close()
        return df

    def stage_upsert(
        self,
        data: pd.DataFrame,
        table: str,
        index_cols: List[str] = None,
        stage_format = "parquet",
    ) -> List[tuple]:
        """
        Stage file and perform an upsert operation.
        :param data: pd.DataFrame
                            DataFame containing data that needs to be uploaded to Snowflake.
        :param table: String
                        Name of table that has to be updated
        :param index_cols: List[str]
                            List of index/primary key columns
        :param stage_format: String
                                Format of the intermediate file for stage/upload.
                                Supported: "parquet" (default) or "csv" (with headers).
        :return: [(<No. of rows inserted>, <<No. of rows updated>)]
        """
        # write data to a parquet file
        temp_file_obj = tempfile.NamedTemporaryFile(suffix='_sf_temp_file')
        temp_fpath = temp_file_obj.name

        if stage_format == "parquet":
            data.to_parquet(temp_fpath)
        elif stage_format == "csv":
            # Export to csv, no index
            data.to_csv(temp_fpath, index=False)
        else:
            raise ValueError("stage_format = {} not supported.".format(stage_format))

        # col_names = ParquetFile(temp_fpath).schema.names
        col_names = list(filter(None, data.columns.tolist() + [idx for idx in data.index.names]))

        # Not sure if we should del data to free up memory. This is going to be a very large dataset.
        # del data

        # check if df columns same as table columns
        query = ("DESC TABLE {tb}".format(tb=table))

        cs = self.con.cursor()
        try:
            desc_result = pd.DataFrame(cs.execute(query))

        except Exception:
            logging.exception("Error in query {q} \n".format(q=query))
            return False
        finally:
            cs.close()

        table_cols = list(desc_result[0])
        col_types = list(desc_result[1])
        # 4th field says "N" if column is "NOT NULL" or "Y" otherwise
        col_allow_null = list(desc_result[3])
        col_allow_null_str = ["NOT NULL" if x=="N" else "" for x in col_allow_null]

        if not index_cols:
            index_cols = list(desc_result[desc_result[5] == 'Y'][0])
        else:
            # Check if user specified index columns exist in Snowflake table
            if set(index_cols).issubset(set(table_cols)) is False:
                raise ValueError("index_cols not present in table.")

        # Check if columns names match columns in Snowflake table
        if set(col_names) != set(table_cols):
            raise ValueError("DataFrame schema doesn't match table schema")

        # create or replace temporary stage
        # "or replace" is useful in bulk uploads with one connection, in order to avoid
        # potential "snowflake.connector.errors.ProgrammingError: ... SQL compilation
        # error: Object 'TEMPSTAGE577' already exists.
        stage_name = "tempstage{}".format(str(randint(100, 999)))

        # Use the default PARQUET_FORMAT
        staging_table_query = (
            "create or replace temporary stage {stage_name} \
            file_format=PARQUET_FORMAT;".format(stage_name=stage_name)
        )
        if stage_format == "csv":
            # Assume usual conditions: 1st line header, fields enclosed by double-quote
            staging_table_query = (
                "create or replace temporary stage {stage_name} \
                file_format = ( TYPE = CSV \
                    SKIP_HEADER = 1 FIELD_OPTIONALLY_ENCLOSED_BY = '\"');".format(stage_name=stage_name)
            )

        cs = self.con.cursor()
        try:
            cs.execute(staging_table_query)
        except Exception as e:
            logging.exception("Error in creating staging table: \n {}".format(e))
            return False
        finally:
            cs.close()

        staging_query = "PUT file:///{fp} @{ts};".format(fp=str(temp_fpath), ts=stage_name)
        cs = self.con.cursor()
        try:
            cs.execute(staging_query)
        except Exception as e:
            logging.exception("Error in staging query: \n {}".format(e))
            return False
        finally:
            cs.close()

        # create temporary table
        temp_table = "temptable{}".format(str(randint(100, 999)))
        temp_schema = ", ".join(
            "{} {} {} ".format(col, typ, null_str)
            for col, typ, null_str in zip(table_cols, col_types, col_allow_null_str)
        )
        temp_tbl_query = 'create or replace temporary table ' \
                         '{tmp_tbl} ({schema});'.format(tmp_tbl=temp_table, schema=temp_schema)

        cs = self.con.cursor()
        try:
            cs.execute(temp_tbl_query)
        except Exception as e:
            logging.exception("Error in creating temporary table: \n {}".format(e))
            return False
        finally:
            cs.close()

        # COPY data into temp table
        copy_query = 'COPY INTO {tmp_tbl} FROM ' \
                     '(SELECT {schema} ' \
                     'FROM @{stg}) ;'.format(tmp_tbl=temp_table,
                                             schema=", ".join("$1:{}::{}".format(col, typ)
                                                              for col, typ in zip(table_cols, col_types)),
                                             stg=stage_name)
        if stage_format == "csv":
            copy_query = (
                'COPY INTO {tmp_tbl} FROM @{stg}/{fp};'.format(
                    tmp_tbl=temp_table, stg=stage_name, fp=Path(temp_fpath).name
                )
            )

        cs = self.con.cursor()
        try:
            cs.execute(copy_query)
        except Exception as e:
            logging.exception("Error in copying data to temporary table: \n {}".format(e))
            return False
        finally:
            cs.close()

        non_idx_cols = list(set(table_cols) - set(index_cols))

        # creating join string
        join_string = " AND ".join("{}.{} = {}.{}".format(table, idx, stage_name, idx) for idx in index_cols)

        # creating UPDATE SET string
        update_string = "UPDATE SET {}".format(
            ", ".join('{col} = {stg}.{col} '.format(col=column, stg=stage_name) for column in table_cols))

        insert_string = "INSERT ({}) VALUES ({})".format(", ".join(str(e) for e in table_cols),
                                                         ", ".join(stage_name + "." + str(e) for e in table_cols))

        upsert_query = "MERGE INTO {tbl} USING \
                        (SELECT * \
                        FROM {tmp_tbl}) \
                        {stg} ON {join}\
                        WHEN MATCHED THEN \
                        {updt} \
                        WHEN NOT MATCHED THEN \
                        {insert};".format(tbl=table, tmp_tbl=temp_table, stg=stage_name, join=join_string,
                                          updt=update_string, insert=insert_string)

        cs = self.con.cursor()
        try:
            cs.execute(upsert_query)
            result = cs.fetchall()
        except Exception as e:
            logging.exception("Error in upsert query: \n {}".format(e))
            return False
        finally:
            cs.close()

        del temp_file_obj
        return result


if __name__ == '__main__':
    df = pd.read_csv('/Users/asavla/Documents/REsurety/met_test3.csv')
    CREDENTIALS = get_secret(secret_name='redas/prod', region_name='us-east-1')
    sf = ReSnowflake(warehouse='DATALOADING', database='PRODUCTION_DATAFEEDS', schema='ERA',
                     sf_user=CREDENTIALS['SNOWFLAKE_USER'], sf_password=CREDENTIALS['SNOWFLAKE_PASS'],
                     sf_account=CREDENTIALS['SNOWFLAKE_ACCOUNT'])
    sf.stage_upsert(data=df, table='ERA_RAW') #, index_cols=["TIME", "LATITUDE", "LONGITUDE"])
