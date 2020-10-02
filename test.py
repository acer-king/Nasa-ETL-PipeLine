import sys
import netCDF4 as nc4
import numpy as np
import urllib3
from typing import List
import time
import cdsapi
import logging
import pandas as pd
import tempfile
from merra2.nasaingestion import get_era5_netcdf
import pandas as pd
import os
import datetime

df = get_era5_netcdf(year=2020, month=8, geo_subset=[
                     35, -106, 29, -98], era5_variables=['v10', 't2m'])
# # del df['t2m']
df = df[:3]
print(df)

# df = get_era5_netcdf(year=2020, month=8, geo_subset=[
#                      35, -106, 29, -98], era5_variables=['v10', 't2m'])
# # del df['t2m']
# df = df[:3]
# print(df)
# 532950
