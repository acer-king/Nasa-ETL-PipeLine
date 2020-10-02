import unittest
import datetime
from merra2.nasaingestion import get_era5_netcdf
from dateutil.relativedelta import relativedelta


class TestEra5(unittest.TestCase):

    def setUp(self) -> None:
        self.test_year = datetime.datetime.now().year
        self.test_month = 8
        self.test_geo_subset = [35, -106, 29, -98]
        self.test_era5_variables = [
            't2m',
            'v10',
            'u10',
            'v100',
            'u100'
        ]
        self.test_invalid_era5_variables = [
            '5tm',
        ]
        self.test_CDS_API_KEY = '61936:f89a7bdd-9e05-4595-872a-8c09077255c1'
        self.test_Invalid_CDS_API_KEY = '61936:f89a7bdd-9e05-4595-872a-8c09077235c1'

    def test_era5_ingestion(self):
        with self.assertRaises(Exception):
            get_era5_netcdf(self.test_year+1,
                            self.test_month,
                            self.test_geo_subset,
                            self.test_era5_variables,
                            self.test_CDS_API_KEY)
        with self.assertRaises(Exception):
            get_era5_netcdf(self.test_year+1,
                            self.test_month,
                            self.test_geo_subset,
                            self.test_era5_variables,
                            self.test_Invalid_CDS_API_KEY)
        with self.assertRaises(Exception):
            get_era5_netcdf(self.test_year,
                            self.test_month,
                            self.test_geo_subset,
                            self.test_invalid_era5_variables,
                            self.test_Invalid_CDS_API_KEY)

    def test_Result(self):
        df = get_era5_netcdf(self.test_year,
                             self.test_month,
                             self.test_geo_subset,
                             self.test_era5_variables,
                             self.test_CDS_API_KEY)
        self.assertEqual(len(df), 613800)
