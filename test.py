import cdsapi

# UID
# 61936
# API Key
# f89a7bdd-9e05-4595-872a-8c09077255c1
x = (0,)
print(len(x))
exit(0)
cds = cdsapi.Client(key="f89a7bdd-9e05-4595-872a-8c09077255c1", url="https://cds.climate.copernicus.eu/api/v2")
cds.retrieve('reanalysis-era5-pressure-levels', {
           "variable": "temperature",
           "pressure_level": "1000",
           "product_type": "reanalysis",
           "date": "2017-12-01/2017-12-31",
           "time": "12:00",
           "format": "grib"
       }, 'download.grib')
