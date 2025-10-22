import xarray as xr
import numpy as np
import glob
import matplotlib.pyplot as plt
import pandas as pd
import ordpy as od


#to plot maps
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

#select variable to analyze
var = 'CHL' #variable to analyze
var = 'RRS412' #variable to analyze

#read a netcdf file to get the lat and lon values
data_files = glob.glob('/g100_scratch/userexternal/gocchipi/GLOBAL/tmp_daily/cmems_????-*.nc')

ds = xr.open_dataset(data_files[0],decode_times=False)

#restrict latitutudes as the regional subset, e.g. between 10 and 20 degrees and longitudes between -15 -60 degrees
ds = ds.sel(latitude=slice(10, 20), longitude=slice(-60, -15))
lon = ds.longitude.values[:]
lat = ds.latitude.values[:]
lon2d, lat2d = np.meshgrid(lon, lat)


#read csv file with lyapunov, PE, C

if var == 'CHL':
    data = pd.read_csv('REGIONAL/merged.csv', header=None)
elif var == 'RRS412':
    data = pd.read_csv('REGIONAL/RRS412_merged.csv', header=None)

data = data.to_numpy()

data[data==0] = np.nan


#the columns of the csv file are: {lat},{lon},{lyap_3_10},{lyap_3_20},{lyap_3_30},{lyap_6_10},{lyap_6_20},{lyap_6_30},{lyap_10_10},{lyap_10_20},{lyap_10_30}
#{coord[0]},{coord[1]},{lyap1},{lyap10},{lyap30},{lyapnolds1},{lyapnolds10},{lyapnolds30},{pe1},{c1},{pe10},{c10},{pe30},{c30}
lat_csv = data[:, 0]
lon_csv = data[:, 1]
lyap_3_10_csv = data[:, 2]
lyap_3_20_csv = data[:, 3]
lyap_3_30_csv = data[:, 4]
lyap_6_10_csv = data[:, 5]
lyap_6_20_csv = data[:, 6]
lyap_6_30_csv = data[:, 7]
lyap_10_10_csv = data[:, 8]
lyap_10_20_csv = data[:, 10]
lyap_10_30_csv = data[:, 12]
lyap_3_180_csv = data[:, 14]
lyap_3_360_csv = data[:, 15]
nolds_3_180_csv = data[:, 16]
nolds_3_360_csv = data[:, 17]
pe_3_180_csv = data[:, 18]
pe_3_360_csv = data[:, 19]
c_3_1_csv = data[:,9]
c_3_10_csv = data[:,11]
c_3_30_csv = data[:,13]
c_3_180_csv = data[:,20]
c_3_360_csv = data[:,21]


# Create 2D grids for lon and lat
lon2d_csv,lat2d_csv = np.meshgrid(lon_csv, lat_csv)
# Create a 2D grid for the lyap, PE and C values
lyap_3_10_2d_csv, _ = np.meshgrid(lyap_3_10_csv, lat_csv)
lyap_3_20_2d_csv, _ = np.meshgrid(lyap_3_20_csv, lat_csv)
lyap_3_30_2d_csv, _ = np.meshgrid(lyap_3_30_csv, lat_csv)
lyap_6_10_2d_csv, _ = np.meshgrid(lyap_6_10_csv, lat_csv)
lyap_6_20_2d_csv, _ = np.meshgrid(lyap_6_20_csv, lat_csv)
lyap_6_30_2d_csv, _ = np.meshgrid(lyap_6_30_csv, lat_csv)
lyap_10_10_2d_csv, _ = np.meshgrid(lyap_10_10_csv, lat_csv)
lyap_10_20_2d_csv, _ = np.meshgrid(lyap_10_20_csv, lat_csv)
lyap_10_30_2d_csv, _ = np.meshgrid(lyap_10_30_csv, lat_csv)
lyap_3_180_2d_csv, _ = np.meshgrid(lyap_3_180_csv, lat_csv)
lyap_3_360_2d_csv, _ = np.meshgrid(lyap_3_360_csv, lat_csv)
nolds_3_180_2d_csv, _ = np.meshgrid(nolds_3_180_csv, lat_csv)
nolds_3_360_2d_csv, _ = np.meshgrid(nolds_3_360_csv, lat_csv)
pe_3_180_2d_csv, _ = np.meshgrid(pe_3_180_csv, lat_csv)
pe_3_360_2d_csv, _ = np.meshgrid(pe_3_360_csv, lat_csv)
c_3_1_2d_csv, _ = np.meshgrid(c_3_1_csv, lat_csv)
c_3_10_2d_csv, _ = np.meshgrid(c_3_10_csv, lat_csv)
c_3_30_2d_csv, _ = np.meshgrid(c_3_30_csv, lat_csv)
c_3_180_2d_csv, _ = np.meshgrid(c_3_180_csv, lat_csv)
c_3_360_2d_csv, _ = np.meshgrid(c_3_360_csv, lat_csv)


hmax, cmax = od.maximum_complexity_entropy(dx=6).T
hmin, cmin = od.minimum_complexity_entropy(dx=6, size=719).T



#add subplots in a grid
fig, axs = plt.subplots(5, 5, figsize=(20, 16), subplot_kw={'projection': ccrs.PlateCarree()})
# Loop through each subplot and plot the data
for i, ax in enumerate(axs.flat):
    ax.coastlines()
    masked = np.zeros_like(lon2d)
    if i == 0:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_3_10_2d_csv, cmap='viridis')
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_3_10_csv, cmap='viridis', s=4)
        ax.set_title('Lyapunov Exponent dim=3 tau=1 day')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 1:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_3_20_2d_csv, cmap='viridis')
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_3_20_csv, cmap='viridis', s=4)
        ax.set_title('Lyapunov Exponent dim=3 tau=10 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 2:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_3_30_2d_csv, cmap='viridis')
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_3_30_csv, cmap='viridis', s=4)
        ax.set_title('Lyapunov Exponent dim=3 tau=30 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 3:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_3_180_2d_csv, cmap='viridis')
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_3_180_csv, cmap='viridis', s=4)
        ax.set_title('Lyapunov Exponent dim=3 tau=180 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 4:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_3_360_2d_csv, cmap='viridis')
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_3_360_csv, cmap='viridis', s=4)
        ax.set_title('Lyapunov Exponent dim=3 tau=360 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 5:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_6_10_2d_csv, cmap='viridis')
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_6_10_csv, cmap='viridis', s=4)
        ax.set_title('Lyapunov NOLDS dim=3 tau=1 day')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 6:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_6_20_2d_csv, cmap='viridis')
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_6_20_csv, cmap='viridis', s=4)
        ax.set_title('Lyapunov NOLDS dim=3 tau=10 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 7:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_6_30_2d_csv, cmap='viridis')
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_6_30_csv, cmap='viridis', s=4)
        ax.set_title('Lyapunov NOLDS dim=3 tau=30 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 8:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, nolds_3_180_2d_csv, cmap='viridis')
        mm = ax.scatter(lon_csv, lat_csv, c=nolds_3_180_csv, cmap='viridis', s=4)
        ax.set_title('Lyapunov NOLDS dim=3 tau=180 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 9:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, nolds_3_360_2d_csv, cmap='viridis')
        mm = ax.scatter(lon_csv, lat_csv, c=nolds_3_360_csv, cmap='viridis', s=4)
        ax.set_title('Lyapunov NOLDS dim=3 tau=360 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 10:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_10_10_2d_csv, cmap='viridis', vmin=0.5, vmax=1)
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_10_10_csv, cmap='viridis', s=4)
        ax.set_title('Entropy tau=1 day')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 11:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_10_20_2d_csv, cmap='viridis', vmin=0.5, vmax=1)
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_10_20_csv, cmap='viridis',  s=4)
        ax.set_title('Entropy tau=10 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 12:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, lyap_10_30_2d_csv, cmap='viridis', vmin=0.5, vmax=1)
        mm = ax.scatter(lon_csv, lat_csv, c=lyap_10_30_csv, cmap='viridis',  s=4)
        ax.set_title('Entropy tau=30 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 13:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, pe_3_180_2d_csv, cmap='viridis', vmin=0.5, vmax=1)
        mm = ax.scatter(lon_csv, lat_csv, c=pe_3_180_csv, cmap='viridis',  s=4)
        ax.set_title('Entropy tau=180 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 14:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, pe_3_360_2d_csv, cmap='viridis', vmin=0.5, vmax=1)
        mm = ax.scatter(lon_csv, lat_csv, c=pe_3_360_csv, cmap='viridis',  s=4)
        ax.set_title('Entropy tau=360 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 15:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, c_3_1_2d_csv, cmap='viridis', vmin=0., vmax=0.5)
        mm = ax.scatter(lon_csv, lat_csv, c=c_3_1_csv, cmap='viridis',  s=4)
        ax.set_title('Complexity tau=1 day')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 16:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, c_3_10_2d_csv, cmap='viridis', vmin=0., vmax=0.5)
        mm = ax.scatter(lon_csv, lat_csv, c=c_3_10_csv, cmap='viridis',  s=4)
        ax.set_title('Complexity tau=10 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 17:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, c_3_30_2d_csv, cmap='viridis', vmin=0., vmax=0.5)
        mm = ax.scatter(lon_csv, lat_csv, c=c_3_30_csv, cmap='viridis',  s=4)
        ax.set_title('Complexity tau=30 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 18:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, c_3_180_2d_csv, cmap='viridis', vmin=0., vmax=0.5)
        mm = ax.scatter(lon_csv, lat_csv, c=c_3_180_csv, cmap='viridis',  s=4)
        ax.set_title('Complexity tau=180 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    elif i == 19:
        #mm = ax.pcolormesh(lon2d_csv, lat2d_csv, c_3_360_2d_csv, cmap='viridis', vmin=0., vmax=0.5)
        mm = ax.scatter(lon_csv, lat_csv, c=c_3_360_csv, cmap='viridis',  s=4)
        ax.set_title('Complexity tau=360 days')
        cbar = plt.colorbar(mm, ax=ax, orientation='vertical', shrink=0.4)
    #plot the complexity-entropy plane in the last 5 subplots
    elif i == 20:
        ax.plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
        ax.plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
        ax.scatter(lyap_10_10_csv, c_3_1_csv, s=4, c='k')
        ax.set_xlabel('Entropy')
        ax.set_ylabel('Complexity')
        ax.set_title('tau = 1 day')
    elif i == 21:
        ax.plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
        ax.plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
        ax.scatter(lyap_10_20_csv, c_3_10_csv, s=4, c='k')
        ax.set_xlabel('Entropy')
        ax.set_ylabel('Complexity')
        ax.set_title('tau = 10 days')
    elif i == 22:
        ax.plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
        ax.plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
        ax.scatter(lyap_10_30_csv, c_3_30_csv, s=4, c='k')
        ax.set_xlabel('Entropy')
        ax.set_ylabel('Complexity')
        ax.set_title('tau = 30 days')
    elif i == 23:
        ax.plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
        ax.plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
        ax.scatter(pe_3_180_csv, c_3_180_csv, s=4, c='k')
        ax.set_xlabel('Entropy')
        ax.set_ylabel('Complexity')
        ax.set_title('tau = 180 days')
    elif i == 24:
        ax.plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
        ax.plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
        ax.scatter(pe_3_360_csv, c_3_360_csv, s=4, c='k')
        ax.set_xlabel('Entropy')
        ax.set_ylabel('Complexity')
        ax.set_title('tau = 360 days')

    if i < 20:
        # Set the extent for the map
        ax.set_extent([-18, -15, 16, 18], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.LAND, color='lightgray')
        ax.gridlines(draw_labels=False, linestyle='--', alpha=0.1)
        #add major lat and lon ticks
        ax.set_xticks(np.arange(-18, -14, 1), crs=ccrs.PlateCarree())
        ax.set_yticks(np.arange(16, 19, 1), crs=ccrs.PlateCarree())
        ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda x, pos: f'{x:.0f}°W'))
        ax.yaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda y, pos: f'{y:.0f}°N'))
        #color the background of each subplot in black
        ax.set_facecolor('black')

# Add a single colorbar for all subplots
#cbar = plt.colorbar(mm, ax=axs, orientation='vertical', shrink=0.4)
fig.tight_layout()
fig.savefig(f'{var}_regional_test_lyapunov.png', dpi=300)

