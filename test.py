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
from merra2.nasaingestion import *
import pandas as pd
import os


df = get_era5_netcdf()
print(df)
