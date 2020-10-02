import os
import urllib3
from typing import List
import time
import cdsapi
import logging
import pandas as pd
import tempfile
from merra2 import FILE_NAME_ERA
import datetime
import netCDF4 as nc4

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
logging.basicConfig(
    format='%(asctime)s  %(levelname)s: %(message)s', level=logging.INFO)
g_era5_variables = ['2m_temperature',
                    '10m_v_component_of_wind',
                    '10m_u_component_of_wind',
                    '100m_v_component_of_wind',
                    '100m_u_component_of_wind'
                    ]


def fill_datetime_components(dt_component_list, fill_positions, suffix=''):
    """
    Fill datetime components (months, days, hours) of a numeric list to put them in a valid format for the CDS
    API request.
    :param dt_component_list: List containing the datetime component to be treated.
    :param fill_positions: Positions to fill with zeros (numeric).
    :param suffix: Optional string to append at the end of each list element. Useful to add ':00' to numeric hours.
    :return: New list containing the elements of the original list treated as requested.
    """
    if not isinstance(dt_component_list, list):
        dt_component_list = [dt_component_list]

    return [str(dtComponent).zfill(fill_positions) + suffix
            for dtComponent in dt_component_list]


def get_era5_netcdf(year: int = None, month: int = None,
                    geo_subset: list = None,
                    era5_variables: list = None,
                    CDS_API_KEY: str = None,
                    outputPath: str = None) -> str:
    """
    Download a NetCDF file from ERA5 and convert the file to pandas object to be returned. The Climate Data
    Store (CDS) API is required. Setup instructions in
    https://cds.climate.copernicus.eu/api-how-to
    1. Create an account at the CDS registration page.
    2. Install the CDS API key in the file $HOME/.cdsapirc
    3. Install the CDS API client: pip install cdsapi
    4. On the CDS site, agree to the Terms of Use of ERA5 datasets.
    Parameters
    ----------
    :param year: int, default 2020
                    year (numeric) in format YYYY.
    :param month: int, default 8
                    month (numeric, valid values from 1 to 12).
    :param era5_variables: strlist, default ['v10', 'u10', 'v100', 'u100', 't2m']
                    list of variables
    :param geo_subset: list of integer, default [35, -106, 29, -98]
                        (35,-106): top left (lat,longitude)
                        (29,-98): right bottom (lat,longitude)
    :param CDS_API_KEY: String, default '61936:f89a7bdd-9e05-4595-872a-8c09077255c1'
                        CDS_API_KEY.
    :param outputPath: String or None, default None
                        root path  to output of zipped csv file
                        if this param is None, then not output csv file
                        if this param is '', then output zipped csv file current path
    Returns
    ------
    :return: None or Pandas Object
                        Assemble all variables as columns into a single PANDAS table with rows by lat, lon,&
                        timestamp UTC time index
    Examples
    -------
    >>> df =get_era5_netcdf(year=2020, month=8, geo_subset=[35, -106, 29, -98], era5_variables=['v10','t2m'])
    >>> df = df[:3]
    >>> df
    ...
    1057008 35.0 -106.00 -1.642295    302.547598
    1057008 35.0 -105.75 -2.697767    302.859733
    1057008 35.0 -105.50 -3.321036    302.855825
    """

    year = 2020 if year is None else year
    month = 8 if month is None else month
    era5_variables = ['v10', 'u10', 'v100', 'u100',
                      't2m'] if era5_variables is None else era5_variables
    global g_era5_variables
    geo_subset = [35, -106, 29, -98] if geo_subset is None else geo_subset
    CDS_API_KEY = '61936:f89a7bdd-9e05-4595-872a-8c09077255c1' if CDS_API_KEY is None else CDS_API_KEY

    if not isinstance(year, int):
        year = int(year)
    if year < 1950:
        raise Exception('Only year greater or equal than 1950 are allowed.')
    if year > datetime.datetime.now().year:
        raise Exception('year cannot be greater than current year: {}'.format(
            datetime.datetime.now().year))

    if month < 1 or month > 12:
        raise Exception('month has to be between 1 and 12.')

    # Determine bounding boxes
    if geo_subset[0] < geo_subset[2] or geo_subset[0] > 50 or geo_subset[2] < 25 \
            or geo_subset[1] > geo_subset[3] or geo_subset[1] < -125 or geo_subset[3] > -67:
        raise Exception(
            'Invalid geoSubset, Please choose a region within us-continental')
    era5_variables = list(
        set(['v10', 'u10', 'v100', 'u100', 't2m']).intersection(era5_variables))
    if len(era5_variables) == 0:
        raise Exception(
            "Invalid era5_variables. should contain vars in ['v10', 'u10', 'v100', 'u100', 't2m']")

    year_filter = str(year)
    month_filter = fill_datetime_components(month, fill_positions=2)

    # List of all days in a month
    day_list = fill_datetime_components(list(range(1, 32)), fill_positions=2)
    day_filter = fill_datetime_components(day_list, fill_positions=2)

    # List of all the hours in a day
    hour_list = fill_datetime_components(
        list(range(0, 24)), fill_positions=2, suffix=":00")
    hour_filter = fill_datetime_components(
        hour_list, fill_positions=2, suffix=':00')

    # Nasa's raw data file name
    temp_fpath = '_era5_{}_.nc'.format(str(year) + str(month))

    # find out if same name file exists, then load it. if not, download dat to a new file with temp_fpath
    global FILE_NAME_ERA
    if FILE_NAME_ERA == temp_fpath or os.path.exists(temp_fpath):
        df = convertCDN2Panda(temp_fpath)
        for field in ['v10', 'u10', 'v100', 'u100', 't2m']:
            if field in era5_variables:
                pass
            else:
                del df[field]
        compression_opts = dict(method='zip', archive_name='out.csv')
        if outputPath is None:
            return df
        try:
            df.to_csv(os.path.join(outputPath, 'out.zip'),
                      compression=compression_opts)
        except OSError:
            logging.info("can't acess file system")
        finally:
            pass
        return df
    else:
        if os.path.exists(str(FILE_NAME_ERA)):
            try:
                os.remove(temp_fpath)
            except OSError:
                logging.info("The old file could not be deleted")
            finally:
                pass

    FILE_NAME_ERA = temp_fpath
    c = cdsapi.Client(
        key=CDS_API_KEY, url="https://cds.climate.copernicus.eu/api/v2")

    # call nasa's API to get raw data file from
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'variable': g_era5_variables,
            'product_type': 'reanalysis',
            'year': year_filter,
            'month': month_filter,
            'day': day_filter,
            'area': geo_subset,
            'time': hour_filter,
            'format': 'netcdf'
        },
        temp_fpath
    )

    logging.info("ERA5 raw file downloaded for Year: {}, Month: {}".format(
        year, month))

    del c
    df = convertCDN2Panda(temp_fpath)
    for field in ['v10', 'u10', 'v100', 'u100', 't2m']:
        if field in era5_variables:
            pass
        else:
            del df[field]

    compression_opts = dict(method='zip', archive_name='out.csv')
    if(outputPath is None):
        return df
    try:
        df.to_csv(os.path.join(outputPath, 'out.zip'),
                  compression=compression_opts)
    except OSError:
        raise Exception("can't create file")
    finally:
        pass
    return df


def convertCDN2Panda(file_name: str = None) -> object:
    """
    imports a Nasa file and exports contents as pandas object
    Parameters
    ----------
    file_name : str or None, default None
                Path to a file to be imported
    Returns
    -------
    Pandas DataFrame
    """
    if(file_name is None):
        return None
    # read in the data
    merraData = nc4.Dataset(file_name, 'r')
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
        return

    data = {}
    for SDS_NAME in fileVars:
        try:
            # read merra data as a vector
            data[SDS_NAME] = merraData.variables[SDS_NAME][:]

        except Exception:
            print(
                'There is an issue with your MERRA file (might be the wrong MERRA file type). Skipping...')
            continue

    timelist = list(merraData.variables['time'][:])
    # convert timestamp to isotime
    timelist = [str(datetime.datetime.fromtimestamp(
        timestamp).isoformat()) for timestamp in timelist]
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
    return df
