import merra2.nasaingestion as merra
df = merra.convertCDN2Panda('_era5_2020-09-28.nc')
print(df)
