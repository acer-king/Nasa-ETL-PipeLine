#!/usr/bin/python
import datetime
import pyutils

CREDENTIALS = pyutils.get_secret(secret_name='redas/prod', region_name='us-east-1')

S3_BUCKET = 'resurety-nonproprietary-datafeeds'

BOUNDING_BOX = {'us_continental': [50, -125, 25, -67],  # North, West, South, East
                'au_continental': [-10, 112, -42, 160],
                'global': [90, -180, -90, 180],
                'testdata': [34.5, -101.75, 33.75, -101],
                }

ERA5_VARIABLES = {
                # Basic met data:
                # Variable list:
                # https://confluence.ecmwf.int/display/CKB/ERA5%3A+data+documentation#ERA5:datadocumentation-Parameterlistings
                # Note: Variables with 10m wind speed/direction are not useful. See:
                # https://confluence.ecmwf.int/display/CKB/ERA5%3A+daily+ocean+waves+10m+wind+speed+and+direction+renamed+ocean+surface+stress+equivalent+10m+neutral+wind+speed+and+direction
                '100m_u_component_of_wind': 'u100',
                '100m_v_component_of_wind': 'v100',
                '10m_u_component_of_wind': 'u10',
                '10m_v_component_of_wind': 'v10',
                '2m_temperature': 't2m',
                # Solar:
                # https://confluence.ecmwf.int/pages/viewpage.action?pageId=111155332
                # https://confluence.ecmwf.int/download/attachments/52462065/radiation_in_mars.pdf?version=1&modificationDate=1447861155677&api=v2
                # Note: Some fields for solar data are missing in older files
                # For example, (1979-01-01, hours 00 to 06)
                #
                # Note: Surface Radiation Parameters: Joule/m2 to Watt/m2:
                # https://confluence.ecmwf.int/pages/viewpage.action?pageId=104241513
                'surface_net_solar_radiation': 'ssr',
                'surface_solar_radiation_downward_clear_sky': 'ssrdc',
                'surface_solar_radiation_downwards': 'ssrd',
                'total_sky_direct_solar_radiation_at_surface': 'fdir',
                'clear_sky_direct_solar_radiation_at_surface': 'cdir',
                # Additional met data:
                'skin_temperature': 'skt',
                'mean_total_precipitation_rate': 'mtpr',
                'precipitation_type': 'ptype',
                'total_precipitation': 'tp',
                'total_cloud_cover': 'tcc',
                'snow_albedo': 'asn',
                'snow_depth': 'sd',
                'snowfall': 'sf',
                'total_column_water_vapour': 'tcwv',
}

ERA5_RAW_MANIFEST = 'ERA5_RAW_MANIFEST'
ERA5_RAW_TS_TABLE = 'ERA5_RAW'
ERA5_DOWNSCALED_MANIFEST = 'ERA5_WIND_DOWNSCALED_MANIFEST'
ERA5_DOWNSCALED_TS_TABLE = 'ERA5_WIND_DOWNSCALED'


def create_era5_raw_s3_key(year: int, month: int, era5_variable: str, geo_subset: str):
    file_name = "{0:%Y%m}.nc".format(datetime.date(year, month, 1))
    return 'reanalysis/era5/raw/{geo}/{var}/{fname}'.format(geo=geo_subset, var=era5_variable, fname=file_name)


def get_era5_snowflake_conn():
    return pyutils.re_snowflake.ReSnowflake(
        warehouse="DATALOADING",
        database="PRODUCTION_DATAFEEDS",
        role="READ_WRITE_RL",
        schema="ERA",
        sf_user=CREDENTIALS['SNOWFLAKE_USER'],
        sf_password=CREDENTIALS['SNOWFLAKE_PASS'],
        sf_account=CREDENTIALS['SNOWFLAKE_ACCOUNT'],
    )
