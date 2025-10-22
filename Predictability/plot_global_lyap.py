import xarray as xr
import numpy as np
import lyapunovV
import glob
import matplotlib.pyplot as plt
import ordpy as od


#to plot maps
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#choose variable to process
var = 'CHL' #variable to analyze
#var = 'RRS412' #variable to analyze

#read data files to get lon lat
data_files = glob.glob('/g100_scratch/userexternal/gocchipi/GLOBAL/tmp_daily/cmems_????-*.nc')
ds = xr.open_dataset(data_files[0],decode_times=False)
lon = ds.longitude.values[:]
lat = ds.latitude.values[:]
lon2d, lat2d = np.meshgrid(lon, lat)

#read csv file with lyapunov, PE, C
data = np.loadtxt('CSV/merged_lyapunov.csv', delimiter=',')
#remove repeated data where a is the same
unique_indices = np.unique(data[:, 0] + data[:, 1], return_index=True)[1]
data = data[unique_indices]
lat_csv = data[:, 0]
lon_csv = data[:, 1]
lyap_csv = data[:, 2]
PE_csv = data[:, 3]
C_csv = data[:, 4]


#lon = ds.lon.values[:]
#lat = ds.lat.values[:]
#lon2d, lat2d = np.meshgrid(lon, lat)

#plot lyapunov map

plt.figure(figsize=(13,6.2))
axins = plt.subplot(111, projection=ccrs.PlateCarree())
axins.coastlines()
masked = np.zeros_like(lon2d)
    # Plot the black overlay in the inset
axins.pcolormesh(lon2d, lat2d, masked, cmap='Reds',vmin=0, vmax=1)#, transform=ccrs.PlateCarree())
axins.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
 # Optional: Add other layers or gridlines
axins.add_feature(cfeature.LAND, color='lightgray')
axins.gridlines(draw_labels=False, linestyle='--', alpha=0.5)
#scatter plot of data from csv file
axins.scatter(lon_csv, lat_csv, c=lyap_csv, cmap='viridis', s=10, transform=ccrs.PlateCarree(), label='Lyapunov Exponent')
# Add colorbar
cbar = plt.colorbar(axins.collections[0], ax=axins, orientation='vertical', shrink=0.6)
cbar.set_label('Lyapunov Exponent', fontsize=12)

plt.savefig('lyapunov_exponents_map.png', dpi=300)

# plot PE map

plt.figure(figsize=(13,6.2))
axins = plt.subplot(111, projection=ccrs.PlateCarree())
axins.coastlines()
masked = np.zeros_like(lon2d)
    # Plot the black overlay in the inset
axins.pcolormesh(lon2d, lat2d, masked, cmap='Reds',vmin=0, vmax=1)#, transform=ccrs.PlateCarree())
axins.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
 # Optional: Add other layers or gridlines
axins.add_feature(cfeature.LAND, color='lightgray')
axins.gridlines(draw_labels=False, linestyle='--', alpha=0.5)
#scatter plot of data from csv file
axins.scatter(lon_csv, lat_csv, c=PE_csv, cmap='viridis', s=2, transform=ccrs.PlateCarree(), label='Lyapunov Exponent')
# Add colorbar
cbar = plt.colorbar(axins.collections[0], ax=axins, orientation='vertical', shrink=0.6)
cbar.set_label('PE', fontsize=12)

plt.savefig('PE_exponents_map.png', dpi=300)

# plot C map

plt.figure(figsize=(13,6.2))
axins = plt.subplot(111, projection=ccrs.PlateCarree())
axins.coastlines()
masked = np.zeros_like(lon2d)
    # Plot the black overlay in the inset
axins.pcolormesh(lon2d, lat2d, masked, cmap='Reds',vmin=0, vmax=1)#, transform=ccrs.PlateCarree())
axins.set_extent([-180, 180, -90, 90], crs=ccrs.PlateCarree())
 # Optional: Add other layers or gridlines
axins.add_feature(cfeature.LAND, color='lightgray')
axins.gridlines(draw_labels=False, linestyle='--', alpha=0.5)
#scatter plot of data from csv file
axins.scatter(lon_csv, lat_csv, c=C_csv, cmap='viridis', s=10, transform=ccrs.PlateCarree(), label='Lyapunov Exponent')
# Add colorbar
cbar = plt.colorbar(axins.collections[0], ax=axins, orientation='vertical', shrink=0.6)
cbar.set_label('C', fontsize=12)

plt.savefig('C_exponents_map.png', dpi=300)


