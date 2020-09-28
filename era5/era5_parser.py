from era5 import *
from pyutils import s3_object_exists
import argparse
import boto3
import datetime
from decimal import Decimal
from functools import reduce
import logging
import numpy as np
import netCDF4
import pandas as pd
import tempfile
import time


logging.basicConfig(format='%(asctime)s  %(levelname)s: %(message)s', level=logging.INFO)


def cartesian_product(*arrays):
    """
    Cartesian product of a variable number of arrays.
    :param arrays: Variable number of arrays
    :return: List[tuples] which is a cartesian product of the input arrays.
    """
    la = len(arrays)
    dtype = np.result_type(*arrays)
    arr = np.empty([len(a) for a in arrays] + [la], dtype=dtype)
    for i, a in enumerate(np.ix_(*arrays)):
        arr[..., i] = a
    return arr.reshape(-1, la)


def extract_netcdf(year: int, month: int, geo_subset: str) -> bool:
    """
    Parse ERA5 raw netcdf files
    :param year: int
                    year (numeric) in format YYYY.
    :param month: int
                    month (numeric, valid values from 1 to 12).

    :param geo_subset: Type of geo subset. Currently 'us_continental' (Default), 'au_continental', 'global' and
                        'testdata' are supported
    :return: Boolean
                True: if parsing and upload succeeds.
                False: if process fails.
    """
    time_era = None
    met_list_df = []
    sf = get_era5_snowflake_conn()
    era_raw_manifest = sf.read_sql(query="SELECT  LATITUDE, LONGITUDE \
                                         FROM {table} WHERE GEO_SUBSET='{geo}'".format(table=ERA5_RAW_MANIFEST,
                                                                                       geo=geo_subset))

    lat_long_pairs = [(row['LATITUDE'], row['LONGITUDE']) for row_num, row in era_raw_manifest.iterrows()]

    for long_name, short_name in ERA5_VARIABLES.items():
        s3_key = create_era5_raw_s3_key(year=year, month=month, era5_variable=long_name, geo_subset=geo_subset)
        s3 = boto3.client('s3')
        file_obj = s3.get_object(Bucket=S3_BUCKET, Key=s3_key)
        data = file_obj['Body'].read()
        era = netCDF4.Dataset(s3_key, mode='r', memory=data)

        # Create time index
        if time_era is None:
            ncdf_time_var = era.variables["time"]
            ncdf_time = netCDF4.num2date(ncdf_time_var[:].data, ncdf_time_var.units, only_use_cftime_datetimes=False)
            time_era = pd.to_datetime(pd.Series(ncdf_time), utc=True)

        met_var = []
        for lat, long in lat_long_pairs:
            def get_era_idx(netcdf_dataset, x, y):
                def get_idx(variable, value):
                    idx = np.where(np.array(netcdf_dataset.variables[variable])==value)
                    if len(idx) != 1 or len(idx[0]) != 1:
                        raise ValueError("wrong number of matches")
                    return int(idx[0][0])
                lon_idx = get_idx("longitude", x)
                lat_idx = get_idx("latitude", y)
                return {
                   "lon_idx": lon_idx,
                   "lat_idx": lat_idx
                }
            try:
                lon_lat_idx = get_era_idx(netcdf_dataset = era, x = round(long, 6), y = round(lat, 6))

                temp_met = pd.DataFrame({'LATITUDE': Decimal("{:0.6f}".format(round(lat, 6))),
                                         'LONGITUDE': Decimal("{:0.6f}".format(round(long, 6))),
                                         'TIME': time_era,
                                         short_name.upper(): era[short_name][:, lon_lat_idx["lat_idx"], lon_lat_idx["lon_idx"]]})
            except Exception:
                continue
            met_var.append(temp_met)

        if len(met_var) == 0:
            met_var = [pd.DataFrame(columns=['LATITUDE', 'LONGITUDE', 'TIME', short_name.upper()])]

        merged_df = pd.concat(met_var)
        merged_df.set_index(keys=['LATITUDE', 'LONGITUDE', 'TIME'], drop=True, inplace=True)
        met_list_df.append(merged_df)
        logging.info("ERA5 Parsed for Year: {}, Month: {}, geo_subset: {}, Variable: {}".format(
            year, month, geo_subset, short_name))

    met = reduce(
        lambda x, y: pd.merge(x, y, left_index=True, right_index=True, how='outer', copy=False, validate='1:1'),
        met_list_df)
    return met


def era5_netcdf_to_snowflake(year: int, month: int, geo_subset: str):
    """
    Parse ERA5 raw netcdf files and upload to Snowflake
    :param year: int
                    year (numeric) in format YYYY.
    :param month: int
                    month (numeric, valid values from 1 to 12).

    :param geo_subset: String
                        Type of geo subset. Currently 'us_continental' (Default), 'au_continental', 'global' and
                        'testdata' are supported.
    :return: pd.DataFrame
    """

    start_time = time.time()
    if year < 1950:
        raise Exception('Only year greater or equal than 1950 are allowed.')
    if year > datetime.datetime.now().year:
        raise Exception('year cannot be greater than current year: {}'.format(datetime.datetime.now().year))

    if month < 1 or month > 12:
        raise Exception('month has to be between 1 and 12.')

    # Determine bounding boxes
    if geo_subset not in BOUNDING_BOX.keys():
        raise Exception(
            'Invalid geoSubset: {}, only "us_continental", "au_continental", '
            '"global" and "testdata" supported.'.format(geo_subset))

    # Make sure files for all variables have already been downloaded for requested year, month, geo_subset
    for var in list(ERA5_VARIABLES.keys()):
        s3_key = create_era5_raw_s3_key(year=year, month=month, era5_variable=var, geo_subset=geo_subset)
        missing_variables = []
        if s3_object_exists(bucket=S3_BUCKET, key=s3_key) is False:
            missing_variables.append(var)

        if missing_variables:
            raise Exception("Download raw netCDF file for year: {}, month:{}, geo:subset: {} \n variables: {}".format(
                year, month, geo_subset, missing_variables))

    met = extract_netcdf(year=year, month=month, geo_subset=geo_subset)

    # ordering columns
    met = met[["T2M", "U100", "V100", "U10", "V10", "SSR", "SSRD", "FDIR", "CDIR", "SSRDC", "ASN", "MTPR",
               "PTYPE", "TCC", "TP", "SD", "SF", "SKT", "TCWV"]]

    sf = get_era5_snowflake_conn()
    try:
        sf.stage_upsert(data=met, table=ERA5_RAW_TS_TABLE, index_cols=["TIME", "LATITUDE", "LONGITUDE"])
    except Exception as e:
        logging.error("ERA5 data upload to snowflake failed for year: {}, month: {}, geo_subset: {}. \n {}".format(
            year, month, geo_subset, e))
        return False

    # deleting temporary file

    end_time = time.time()
    logging.info("ERA5 Parsed for Year: {}, Month: {}, geo_subset: {}. \n Time taken: {}".format(
        year, month, geo_subset, end_time - start_time))
    return True


def era5_netcdf_to_s3(year: int, month: int, geo_subset: str):
    """
    Parse ERA5 raw netcdf files and upload to s3
    :param year: int
                    year (numeric) in format YYYY.
    :param month: int
                    month (numeric, valid values from 1 to 12).

    :param geo_subset: String
                        Type of geo subset. Currently 'us_continental' (Default), 'au_continental', 'global' and
                        'testdata' are supported.
    :return: pd.DataFrame
    """

    start_time = time.time()
    if year < 1950:
        raise Exception('Only year greater or equal than 1950 are allowed.')
    if year > datetime.datetime.now().year:
        raise Exception('year cannot be greater than current year: {}'.format(datetime.datetime.now().year))

    if month < 1 or month > 12:
        raise Exception('month has to be between 1 and 12.')

    # Determine bounding boxes
    if geo_subset not in BOUNDING_BOX.keys():
        raise Exception(
            'Invalid geoSubset: {}, only "us_continental", "au_continental", '
            '"global" and "testdata" supported.'.format(geo_subset))

    # Make sure files for all variables have already been downloaded for requested year, month, geo_subset
    for var in list(ERA5_VARIABLES.keys()):
        s3_key = create_era5_raw_s3_key(year=year, month=month, era5_variable=var, geo_subset=geo_subset)
        missing_variables = []
        if s3_object_exists(bucket=S3_BUCKET, key=s3_key) is False:
            missing_variables.append(var)

        if missing_variables:
            raise Exception("Download raw netCDF file for year: {}, month:{}, geo:subset: {} \n variables: {}".format(
                year, month, geo_subset, missing_variables))

    met = extract_netcdf(year=year, month=month, geo_subset=geo_subset)

    # ordering columns
    met = met[["T2M", "U100", "V100", "U10", "V10", "SSR", "SSRD", "FDIR", "CDIR", "SSRDC", "ASN", "MTPR",
               "PTYPE", "TCC", "TP", "SD", "SF", "SKT", "TCWV"]]

    file_name = "{0:%Y%m}.parquet".format(datetime.date(year, month, 1))
    
    temp_file_obj = tempfile.NamedTemporaryFile(suffix=file_name)
    temp_fpath = temp_file_obj.name
    met.to_parquet(temp_fpath)

    s3_key = 'reanalysis/era5/processed/{geo}/{fname}'.format(geo=geo_subset, fname=file_name)
    
    s3 = boto3.client('s3')
    s3.upload_file(Filename=temp_fpath, Bucket=S3_BUCKET, Key=s3_key)
    
    # deleting temporary file
    del temp_file_obj

    end_time = time.time()
    logging.info("ERA5 Parsed and uploaded to S3 for Year: {}, Month: {}, geo_subset: {}. \n Time taken: {}".format(
        year, month, geo_subset, end_time - start_time))
    return True


if __name__ == '__main__':
    my_parser = argparse.ArgumentParser(description='Run ERA5 parser')
    my_parser.add_argument('year', metavar='year', type=int, help='Request year')
    my_parser.add_argument('month', metavar='month', type=int, help='Request month')
    my_parser.add_argument('geo_subset', metavar='geo_subset', type=str, help='Request geo_subset')
    my_parser.add_argument('upload_to', metavar='upload_to', type=str, default="snowflake", choices=['snowflake', 's3'], help='Upload to snowflake or s3')

    args = my_parser.parse_args()

    logging.info("Parsing ERA5 netcdf for Year: {}, "
                 "Month: {}, geo_subset: {}, ".format(args.year, args.month, args.geo_subset))
    if args.upload_to == 'snowflake':
        era5_netcdf_to_snowflake(year=args.year, month=args.month, geo_subset=args.geo_subset)
    elif args.upload_to == 's3':
        era5_netcdf_to_s3(year=args.year, month=args.month, geo_subset=args.geo_subset)
