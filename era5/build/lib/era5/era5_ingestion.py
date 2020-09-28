from era5 import *
from pyutils import s3_object_exists
import argparse
import boto3
import cdsapi
import logging
import pandas as pd
import tempfile
import time
import urllib3
from typing import List

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
logging.basicConfig(format='%(asctime)s  %(levelname)s: %(message)s', level=logging.INFO)


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

    return [str(dtComponent).zfill(fill_positions) + suffix \
            for dtComponent in dt_component_list]


def download_era5_netcdf(year: int, month: int, era5_variable: str, geo_subset: str,
                         force_download: bool = False) -> str:
    """
    Download a NetCDF file from ERA5. The Climate Data
    Store (CDS) API is required. Setup instructions in
    https://cds.climate.copernicus.eu/api-how-to
    1. Create an account at the CDS registration page.
    2. Install the CDS API key in the file $HOME/.cdsapirc
    3. Install the CDS API client: pip install cdsapi
    4. On the CDS site, agree to the Terms of Use of ERA5 datasets.
    :param year: int
                    year (numeric) in format YYYY.
    :param month: int
                    month (numeric, valid values from 1 to 12).
    :param era5_variable: str
                            ERA5 variable name. Has to be in ERA5_VARIABLES
    :param geo_subset: String
                        Type of geo subset. Currently, only support 'us_continental' (Default), au_continental
                        'global' and 'testdata'.
    :param force_download: Boolean
                            If True: will download file even if one exists in S3,
                            If False: will skip download if file exists in S3
    :return: String
                S3 URI
    """
    if not isinstance(year, int):
        year = int(year)

    if year < 1950:
        raise Exception('Only year greater or equal than 1950 are allowed.')
    if year > datetime.datetime.now().year:
        raise Exception('year cannot be greater than current year: {}'.format(datetime.datetime.now().year))

    if month < 1 or month > 12:
        raise Exception('month has to be between 1 and 12.')

    # Determine bounding boxes
    if geo_subset not in BOUNDING_BOX.keys():
        raise Exception(
            'Invalid geoSubset, only "us_continental", "au_continental", "global" and "testdata" supported.')
    area_subset = BOUNDING_BOX[geo_subset]
    year_filter = str(year)
    month_filter = fill_datetime_components(month, fill_positions=2)

    # List of all days in a month
    day_list = fill_datetime_components(list(range(1, 32)), fill_positions=2)
    day_filter = fill_datetime_components(day_list, fill_positions=2)

    # List of all the hours in a day
    hour_list = fill_datetime_components(list(range(0, 24)), fill_positions=2, suffix=":00")
    hour_filter = fill_datetime_components(hour_list, fill_positions=2, suffix=':00')

    s3 = boto3.client('s3')
    s3_key = create_era5_raw_s3_key(year=year, month=month, era5_variable=era5_variable, geo_subset=geo_subset)

    # If force_download is False, only download if file doesn't exist in s3.
    if force_download is False:
        if s3_object_exists(bucket=S3_BUCKET, key=s3_key):
            logging.info("ERA5 raw file already exists for Year: {}, Month: {}, Variable: {}".format(
                year, month, era5_variable))
            return True

    temp_file_obj = tempfile.NamedTemporaryFile(suffix='_era5_{}_{}_.nc'.format(str(year) + str(month), era5_variable))
    temp_fpath = temp_file_obj.name

    era5_client_key = CREDENTIALS['ERA5_CDS_API']
    c = cdsapi.Client(key=era5_client_key, url="https://cds.climate.copernicus.eu/api/v2")

    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'variable': [era5_variable],
            'product_type': 'reanalysis',
            'year': year_filter,
            'month': month_filter,
            'day': day_filter,
            'area': area_subset,
            'time': hour_filter,
            'format': 'netcdf'
        },
        temp_fpath
    )

    s3.upload_file(Filename=temp_fpath, Bucket=S3_BUCKET, Key=s3_key)
    logging.info("ERA5 raw file downloaded for Year: {}, Month: {}, Variable: {}".format(
        year, month, era5_variable))
    del temp_file_obj
    del c
    return "s3://{bucket}/{s3_key}".format(bucket=S3_BUCKET, s3_key=s3_key)


def download_era5_netcdf_all_variables(year: int, month: int, geo_subset: str = 'us_continental',
                                       force_download: bool = False) -> List[str]:
    """
    Download NetCDF files from ERA5. The Climate Data
    Store (CDS) API is required. Setup instructions in
    https://cds.climate.copernicus.eu/api-how-to
    1. Create an account at the CDS registration page.
    2. Install the CDS API key in the file $HOME/.cdsapirc
    3. Install the CDS API client: pip install cdsapi
    4. On the CDS site, agree to the Terms of Use of ERA5 datasets.
    :param year: int
                    year (numeric) in format YYYY.
    :param month: int
                    month (numeric, valid values from 1 to 12).
    :param geo_subset: String
                        Type of geo subset. Currently, only support 'us_continental' (Default), 'au_continental',
                        'global' and 'testdata'.
    :param force_download: Boolean
                            If True: will download file even if one exists in S3,
                            If False: will skip download if file exists in S3
    :return: List[str]:
                    List of s3 URI's for the downloaded netcdf files
    """
    start_time = time.time()
    if not isinstance(year, int):
        year = int(year)
    # Jan 2019: Data available from 1979 onwards. By mid 2019, the historic
    # data for 1950-1978 should be available. Consider this on the validation:
    # https://www.ecmwf.int/en/forecasts/datasets/archive-datasets/reanalysis-datasets/era5
    if year < 1950:
        raise Exception('Only year greater or equal than 1950 are allowed.')
    if year > datetime.datetime.now().year:
        raise Exception('year cannot be greater than current year: {}'.format(datetime.datetime.now().year))

    if month < 1 or month > 12:
        raise Exception('month has to be between 1 and 12.')

    # Determine bounding boxes
    if geo_subset not in BOUNDING_BOX.keys():
        raise Exception(
            'Invalid geoSubset, only "us_continental", "au_continental", "global" and "testdata" supported.')
    result_s3_keys_list = []
    for var in ERA5_VARIABLES:
        result_s3_keys_list.append(download_era5_netcdf(year=year, month=month, era5_variable=var, geo_subset=geo_subset,
                             force_download=force_download))

    end_time = time.time()
    logging.info("ERA5 raw files downloaded for Year: {}, Month: {}. Download Time: {} seconds".format(
        year, month, end_time - start_time))

    return result_s3_keys_list


def missing_netcdf(geo_subset: str, start_date: datetime.date, end_date: datetime.date) -> bool:
    """
    Prints if there any any missing ERA5 raw netcdf files for requested geo_subset for a user specified date range.
    :param geo_subset: ERA5 geo subset. One of BOUNDING_BOX.keys()
    :param start_date: Start date of the time range you want to check
    :param end_date: End date of the time range you want to check
    :return: True if missing else False
    """
    date_set = set(
        [x.date() for x in pd.date_range(start=start_date, end=end_date, freq='MS')])
    missing_flag = False
    client = boto3.client('s3')
    for var in list(ERA5_VARIABLES.keys()):
        objects = client.list_objects(Bucket=S3_BUCKET,
                                      Prefix='reanalysis/era5/raw/{}/{}/'.format(geo_subset, var))
        files = [datetime.datetime.strptime((i['Key'].split('/')[-1].replace('.nc', '01')), "%Y%m%d").date() for i in
                 objects['Contents']]
        missing = sorted(date_set - set(files))
        if len(missing) > 0:
            missing_flag = True
            logging.info("Geo Subset: {}, Variable: {} missing netcdf files \n  {}".format(geo_subset, var, missing))

    del client
    if missing_flag:
        return True
    return False


if __name__ == '__main__':
    my_parser = argparse.ArgumentParser(description='Run ERA5 parser')
    my_parser.add_argument('year', metavar='year', type=int, help='Request year')
    my_parser.add_argument('month', metavar='month', type=int, help='Request month')
    my_parser.add_argument('geo_subset', metavar='year', type=str, help='Request geo_subset')
    my_parser.add_argument('force_download', metavar='force_download', type=bool, help='force_download')

    args = my_parser.parse_args()

    logging.info("Downloading ERA5 netcdf for Year: {}, Month: {}, geo_subset: {}, "
                 "force_download: {}".format(args.year, args.month, args.geo_subset, args.force_download))

    download_era5_netcdf_all_variables(year=args.year,
                                       month=args.month,
                                       geo_subset=args.geo_subset,
                                       force_download=args.force_download)

    logging.info("Finished downloading ERA5 netcdf for Year: {}, Month: {}, geo_subset: {}, "
                 "force_download: {}".format(args.year, args.month, args.geo_subset, args.force_download))
