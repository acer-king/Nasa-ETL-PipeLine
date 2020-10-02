## Install Required Packages
    git clone https://github.com/acer-king/Nasa-ETL-PipeLine
    cd Nasa-ETL-PipeLine
    python setup.py install
    pip install -r requirements.txt
Next, include the SDK in your project  by adding:

    from merra2.nasaingestion import get_era5_netcdf

## Authentication 

    Store (CDS) API is required. Setup instructions in
    https://cds.climate.copernicus.eu/api-how-to
    1. Create an account at the CDS registration page.
    2. Install the CDS API key in the file $HOME/.cdsapirc
    3. Install the CDS API client: pip install cdsapi
    4. On the CDS site, agree to the Terms of Use of ERA5 datasets.

## Examples
    If you don't speciyfy output filePath, then it just return DataFrme.
    If filePath specified, then It outputs ziped csv file in the path

    df = get_era5_netcdf(year=2020, month=8, geo_subset=[
                         35, -106, 29, -98], era5_variables=['v10', 't2m'])
    df = df[:3]
    print(df)
    
