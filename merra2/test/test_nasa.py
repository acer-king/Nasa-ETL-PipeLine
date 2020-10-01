import unittest
import datetime
from merra2.nasaingestion import *
from dateutil.relativedelta import relativedelta


class TestEra5(unittest.TestCase):

    def setUp(self) -> None:
        self.test_year = datetime.datetime.now().year
        self.test_month = 8
        self.test_geo_subset = 'testdata'
        self.test_era5_variables = [
            '2m_temperature',
            '10m_v_component_of_wind',
            '10m_u_component_of_wind',
            '100m_v_component_of_wind',
            '100m_u_component_of_wind'
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

    def test_Result(self):
        df = get_era5_netcdf(self.test_year,
                             self.test_month,
                             self.test_geo_subset,
                             self.test_era5_variables,
                             self.test_CDS_API_KEY)
        self.assertEqual(len(df), 613800)
