from mpi4py import MPI
import numpy as np
import xarray as xr
import pandas as pd
import lyapunovV
#import nolds as NOLDS

from os import cpu_count
from matplotlib import pyplot as plt
import ordpy as od
import os

    

#used to run over coordinates
def mpi_pyramidal_complexity_entropy_map(
    data, epsilon,t_steps,DELTA,delta0
):

    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    size  = comm.Get_size()

    n_lat = data.sizes['latitude']
    n_lon = data.sizes['longitude']
    total = n_lat * n_lon

    # Split work among ranks
    all_indices = np.arange(total, dtype=np.int64)
    index_chunks = np.array_split(all_indices, size)
    my_indices = index_chunks[rank]

    # Preallocate local results (columns: idx, lyap_palladin, etc)
    n_local = len(my_indices)
    local_arr = np.empty((n_local, 3), dtype=float) #size of indicators + 1, now 2 indicators 

    for i, idx in enumerate(my_indices):
        lat = idx // n_lon
        lon = idx % n_lon
        
        #extract coordinate values
        time_dim, lat_dim, lon_dim  = "time", "latitude", "longitude"
        xx = data.isel({ lat_dim: lat, lon_dim: lon}).to_array()

        #remove nans
        xx = xx[0,:]
        print(xx.shape)
        xx = xx[~np.isnan(xx)]
        xx = xx.values

        #compute lyapunov palladin
        #python class
        lp = lyapunovV.LYAP(xx)
        _,lyap_guido,_ =lp.lyap_e_paladin(epsilon=epsilon,t_steps=t_steps,Delta=DELTA,delta0=delta0)
        try:
            lyap_wolf = lp.lyap_e(dt=1,ndim=3,evolve=20,tau=10)
        except:
            lyap_wolf = np.nan
        ##add othe computations
        # lyap = f(xx)


        local_arr[i, 0] = idx        # coord idx
        local_arr[i, 1] = lyap_guido # associate coordinate and lyap
        local_arr[i, 2] = lyap_wolf # associate coordinate and lyap
#       local_arr[i, 3] = lyap       # uncomment to add other indicators

    # Gather arrays at root
    gathered = comm.gather(local_arr, root=0)

    if rank == 0:
        # Stack all results
        all_results = np.vstack(gathered)

        # Allocate flat arrays
        lyap_guido_flat   = np.empty(total, dtype=float)
        lyap_wolf_flat   = np.empty(total, dtype=float)
        #lyap_add_flat   = np.empty(total, dtype=float)

        # Fill arrays using indices
        idxs = all_results[:, 0].astype(int)
        lyap_guido_flat[idxs]   = all_results[:, 1]
        lyap_wolf_flat[idxs]   = all_results[:, 2]
       #lyap_add_flat[idxs]  = all_results[:, "]

        # Reshape into (lat, lon)
        lyap_guido   = lyap_guido_flat.reshape((n_lat, n_lon))
        lyap_wolf   = lyap_wolf_flat.reshape((n_lat, n_lon))
#       lyap_add  = lyap_add_flat.reshape((n_lat, n_lon))

        return xr.Dataset(
            {
                'lyap_guido':           (('latitude','longitude'), lyap_guido),
                'lyap_wolf':           (('latitude','longitude'), lyap_wolf),
              #  'lyap_addy':        (('latitude','longitude'), lyap_add),
            },
            coords={
                'latitude':  data.coords['latitude'],
                'longitude': data.coords['longitude'],
            }
        )
    else:
        return None

comm  = MPI.COMM_WORLD
rank  = comm.Get_rank()
size  = comm.Get_size()

print('starting computation',flush=True)

# Set parameters
epsilon = 0.02
t_steps = 500
DELTA = 0.2
delta0 = 0.02

# Load data regional test dataset

if not os.path.exists('regional_dataset_CHL_10_20_-30_-15.nc'):
    import shutil
    shutil.copy2('/g100_scratch/userexternal/gocchipi/GLOBAL/tmp_daily/regional_dataset_CHL_10_20_-30_-15.nc', os.getcwd())
    print('copied datafile in working directory')
data = xr.open_dataset('regional_dataset_CHL_10_20_-30_-15.nc', decode_times=False)

# select a smaller region for rapid testing
data = data.sel(latitude=slice(17, 18), longitude=slice(-18, -17))

units, reference_date = data.time.attrs['units'].split('since')

#convert time to datetime64
data['time'] = pd.date_range(start=reference_date, periods=data.sizes['time'], freq='d')

#clean data removing negative and too high values
data['CHL'] = data['CHL'].where((data['CHL'] > 0) & (data['CHL'] < 100))

zoom = np.linspace(0.1,10,5)
# Set parameters
epsilon = 0.02
t_steps = 90
DELTA = 0.2
delta0 = 0.02
fig,axs = plt.subplots(1, 6, figsize=(16,4))
axs = axs.ravel()
for i,iz in enumerate(zoom):
    epsilon = epsilon *iz
    DELTA = DELTA *iz
    delta0 = delta0 *iz
    print(f'trying deltq0={delta0}')
    outfile = f'test_palladin_{delta0}.nc'
    if os.path.exists(outfile):
        results =xr.open_dataset(outfile)
    else:
        results = mpi_pyramidal_complexity_entropy_map(
        data, epsilon,t_steps,DELTA,delta0
        )
        results.to_netcdf(outfile, mode='w', format='NETCDF4')
    print('computation done')
    if rank == 0:
        results.lyap_guido.plot(ax=axs[i], cmap='viridis')
        axs[i].set_title(f'Palladin, delta0= {delta0:.3f}')

if rank == 0:
    results.lyap_wolf.plot(ax=axs[-1], cmap='viridis')
    axs[-1].set_title('Wolf')
    fig.savefig('test_indicators.png')
    fig.show()    
