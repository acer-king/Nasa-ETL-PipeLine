import urllib3
from typing import List
import time
import cdsapi
import logging
import pandas as pd
import tempfile
from merra2.nasaingestion import *
import os

filePath = download_era5_netcdf(era5_variables='100m_u_component_of_wind',geo_subset='testdata',CDS_API_KEY='61936:f89a7bdd-9e05-4595-872a-8c09077255c1')
filePath = download_era5_netcdf(era5_variables='100m_u_component_of_wind',geo_subset='testdata',CDS_API_KEY='61936:f89a7bdd-9e05-4595-872a-8c09077255c1')
print(filePath)
# cds = cdsapi.Client(key="61936:f89a7bdd-9e05-4595-872a-8c09077255c1", url="https://cds.climate.copernicus.eu/api/v2")
# cds.retrieve('reanalysis-era5-pressure-levels', {
#            "variable": "temperature",
#            "pressure_level": "1000",  
#            "product_type": "reanalysis",
#            "date": "2017-12-01/2017-12-31",
#            "time": "12:00",
#            "format": "grib"
#        }, 'download.grib')


