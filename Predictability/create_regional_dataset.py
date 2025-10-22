import xarray as xr
import glob

#choose which variable to process
var = 'CHL'
var = 'RRS412'
#assign the correct data files path according to the variable
if var == 'CHL':
    data_files = glob.glob('/g100_scratch/userexternal/gocchipi/GLOBAL/tmp_daily/cmems_????-*.nc')
elif var == 'RRS412':
    data_files = glob.glob('/g100_scratch/userexternal/icunico0/DATA_RIFLETTANZA/DATA_RRS412/cmmems_RRS412_????-*.nc')
data_files.sort()

#regional area limits
lat_min, lat_max = 10, 20
lon_min, lon_max = -30, -15

#concatenate daily files to obtain a single dataset within the regional area
def create_regional_dataset(data_files, lat_min, lat_max, lon_min, lon_max):
    regional_datasets = []
    
    for file in data_files:
        ds = xr.open_dataset(file, decode_times=False)
        # Restrict latitudes and longitudes
        ds = ds.sel(latitude=slice(lat_min, lat_max), longitude=slice(lon_min, lon_max))
        regional_datasets.append(ds)
    
    # Concatenate all datasets along the time dimension
    regional_ds = xr.concat(regional_datasets, dim='time')
    return regional_ds

ds = create_regional_dataset(data_files, lat_min, lat_max, lon_min, lon_max)
# Save the regional dataset to a new NetCDF file
output_file = f'regional_dataset_{var}_{lat_min}_{lat_max}_{lon_min}_{lon_max}.nc'
ds.to_netcdf(output_file, mode='w', format='NETCDF4')