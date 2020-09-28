from setuptools import setup

setup(name='pyutils',
      version='0.1',
      description='Python utilities',
      url='',
      author='Anshal Savla',
      author_email='asavla@resurety.com',
      license='REsurety INC.',
      packages=['pyutils', 'pyutils.test'],
      install_requires=['boto3', 'botocore', 'pandas', 'snowflake-connector-python', ],
      zip_safe=False)
