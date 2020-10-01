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

# file = download_era5_netcdf()

#!/usr/bin/python
'''
Module: read_merra_aerosol_and_dump_ascii.py
==========================================================================================
Disclaimer: The code is for demonstration purposes only. Users are responsible to check for accuracy and revise to fit their objective.

Author: Justin Roberts-Pierel, 2015
Organization: NASA ARSET
Purpose: To save a MERRA netCDF4 file in ASCII format. Saves lat, lon, and other SDS dependent on file.

See the README associated with this module for more information.
==========================================================================================
'''

# import necessary modules

# loops through all files listed in the text file
for FILE_NAME in ['_era5_20208_.nc']:
    # read in the data
    merraData = nc4.Dataset(FILE_NAME, 'r')
    variables = set(merraData.variables)
    desiredVariables = set(
        {'u10', 'u100', 'v10', 'v100', 't2m'}
    )
    desiredVariables = set([x.lower() for x in desiredVariables])
    var1 = variables.intersection(desiredVariables)
    desiredVariables = set([x.upper() for x in desiredVariables])
    var2 = variables.intersection(desiredVariables)
    fileVars = list(var1.union(var2))
    if len(fileVars) == 0:
        print('This file contains none of the selected SDS. Skipping...')
        continue

    data = {}
    for SDS_NAME in fileVars:
        try:
            # read merra data as a vector
            data[SDS_NAME] = merraData.variables[SDS_NAME][:]

        except Exception:
            print(
                'There is an issue with your MERRA file (might be the wrong MERRA file type). Skipping...')
            continue
    # exit(0)

    timelist = list(merraData.variables['time'][:])
    latlist = list(merraData.variables['latitude'][:])
    lnglist = list(merraData.variables['longitude'][:])

    df = pd.DataFrame(columns=['v10', 'u10', 'v100', 'u100', 't2m'])

    df['v10'] = data['v10'].reshape(-1).tolist()
    df['u10'] = data['u10'].reshape(-1).tolist()
    df['u100'] = data['u100'].reshape(-1).tolist()
    df['v100'] = data['v100'].reshape(-1).tolist()
    df['t2m'] = data['t2m'].reshape(-1).tolist()
    df['time'] = sum([[value]*len(latlist)*len(lnglist)
                      for value in timelist], [])
    df['lat'] = sum([[value]*len(lnglist)
                     for value in latlist], [])*len(timelist)
    df['lng'] = lnglist*len(timelist)*len(latlist)
    df.set_index(['time', 'lat', 'lng'], inplace=True)
    df.to_csv('1.csv')

print('\nAll valid files have been saved successfully.')
