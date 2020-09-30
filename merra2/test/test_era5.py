import unittest
import datetime
import merra2
from dateutil.relativedelta import relativedelta

class TestEra5(unittest.TestCase):
    def setUp(self) -> None:
        self.test_year = datetime.datetime.now().year
        self.test_month = (datetime.datetime.now() - relativedelta(months=+3)).month
        self.test_geo_subset = 'testdata'

    def test_era5_ingestion (self):
        now = datetime.datetime.now()
        # testing download_era5_netcdf
        self.assertRaises(Exception,
                          merra2.download_era5_netcdf,
                          self.test_year + 1,
                          self.test_month,
                          next(iter(ERA5_VARIABLES)),
                          self.test_geo_subset)

        self.assertRaises(Exception,
                          era5_ingestion.download_era5_netcdf,
                          self.test_year,
                          (now + relativedelta(months=1)).month,
                          next(iter(ERA5_VARIABLES)),
                          self.test_geo_subset)

        self.assertRaises(Exception,
                          era5_ingestion.download_era5_netcdf,
                          self.test_year,
                          self.test_month,
                          next(iter(ERA5_VARIABLES)),
                          'invalid_geo_subset')

        self.assertRaises(Exception,
                          era5_ingestion.download_era5_netcdf,
                          self.test_year,
                          self.test_month,
                          'invalid_era_variable',
                          self.test_geo_subset)

        s3_uri = era5_ingestion.download_era5_netcdf(year=self.test_year,
                                                     month=self.test_month,
                                                     era5_variable=next(iter(ERA5_VARIABLES)),
                                                     geo_subset=self.test_geo_subset,
                                                     force_download=True)

        s3_key = s3_uri.replace('s3://resurety-nonproprietary-datafeeds/', '')
        self.assertTrue(pyutils.s3_object_exists(bucket=S3_BUCKET, key=s3_key),
                        msg='download_era5_netcdf failed to download data.')

        s3 = boto3.resource('s3')
        s3.Object(S3_BUCKET, s3_key).delete()

        # testing  download_era5_netcdf_all_variables
        self.assertRaises(Exception,
                          era5_ingestion.download_era5_netcdf_all_variables,
                          self.test_year + 1,
                          self.test_month,
                          self.test_geo_subset)

        self.assertRaises(Exception,
                          era5_ingestion.download_era5_netcdf_all_variables,
                          self.test_year,
                          (now + relativedelta(months=1)).month,
                          self.test_geo_subset)

        self.assertRaises(Exception,
                          era5_ingestion.download_era5_netcdf_all_variables,
                          self.test_year,
                          self.test_month,
                          'invalid_geo_subset')

    def test_era5_sf_tables(self):
        era5_tables = [ERA5_RAW_MANIFEST, ERA5_RAW_TS_TABLE, ERA5_DOWNSCALED_MANIFEST, ERA5_DOWNSCALED_TS_TABLE]
        sf = get_era5_snowflake_conn()
        missing_tables = []
        for table in era5_tables:
            query = ("DESC TABLE {tb}".format(tb=table))
            desc_result = sf.read_sql(query=query)
            if isinstance(desc_result, int) and desc_result == 1:
                missing_tables.append(table)

        self.assertEqual(0, len(missing_tables),
                         msg="The following Snowflake tables are missing: {}".format(missing_tables))

    def test_data_completemess_netcdf(self):
        self.assertFalse(era5_ingestion.missing_netcdf(geo_subset=list(BOUNDING_BOX.keys())[0],
                                                       start_date=datetime.date(1979,1,1),
                                                       end_date=datetime.date(self.test_year, self.test_month, 1)),
                         msg="Missing netcdf files for geo_subset: {}".format(list(BOUNDING_BOX.keys())[0]))

        self.assertFalse(era5_ingestion.missing_netcdf(geo_subset=list(BOUNDING_BOX.keys())[1],
                                                       start_date=datetime.date(1979, 1, 1),
                                                       end_date=datetime.date(self.test_year, self.test_month, 1)),
                         msg="Missing netcdf files for geo_subset: {}".format(list(BOUNDING_BOX.keys())[1]))
