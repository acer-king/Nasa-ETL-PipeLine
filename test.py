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

    timelist = list(merraData.variables['time'][:])
    latlist = list(merraData.variables['latitude'][:])
    lnglist = list(merraData.variables['longitude'][:])

    result = pd.DataFrame()
    for SDS_NAME in fileVars:
        try:
            # read merra data as a vector
            data = merraData.variables[SDS_NAME][:]
        except:
            print(
                'There is an issue with your MERRA file (might be the wrong MERRA file type). Skipping...')
            continue

        print(len(data))
        print(SDS_NAME)
        print(type(data))
        inp = input("")

        continue
        if len(data.shape) == 4:
            level = data.shape[1]-1
            data = data[0, level, :, :]
        elif len(data.shape) == 3:
            level = data.shape[0]-1
            data = data[level, :, :]
        data = data.ravel()
        # save variable
        output[:, index] = data
        tempOutput.append(SDS_NAME)
        index += 1

    # # stacks the titles on the data array to be saved to file
    # output = np.row_stack((tempOutput, output))
    # # create the name of the file from the filename
    # outputName = '{0}.txt'.format(FILE_NAME[:-4])
    # # save the file
    # np.savetxt(outputName, output, delimiter=',', fmt='%s')
print('\nAll valid files have been saved successfully.')
