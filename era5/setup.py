from setuptools import setup

setup(name='era5',
      version='0.1',
      description='Data ingestion and parsing for ERA5',
      url='',
      author='Anshal Savla',
      author_email='asavla@resurety.com',
      license='REsurety INC.',
      packages=['era5', 'era5.test'],
      install_requires=['cdsapi', 'netCDF4', 'numpy', 'pandas', ],
      zip_safe=False)
