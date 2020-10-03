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


if __name__ == "__main__":
    df = get_era5_netcdf(year=2020, month=8, geo_subset=[
                         35, -106, 29, -98], outputPath="", forceUpdate=False)

    print(df[:3])
    print(df[-3:])
    pass
