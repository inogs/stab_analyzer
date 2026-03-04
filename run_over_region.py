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

# Description of functions:
# SAMPLED_LAT_LON: extract timeseries from a spatio-temporal random walk starting from given lat, lon, time
# SAMPLED_HC_LAT_LON: calculate permutation entropy, complexity and lyapunov exponent for a given lat, lon, starting time using random-walk sampling
# mpi_pyramidal_complexity_entropy_map: run in parrallel (mpi4py) over the whole dataset SAMPLED_HC_LAT_LON to produce maps of permutation entropy and complexity

def SAMPLED_LAT_LON(data, lat_center, lon_center, initial_time_layer, steps=10, ensemble=10, p_back=0.2, connectivity=10, 
                       avoid_nan=True, d=3, tau=10):
    """
    Extract a timeseries from a spatio-temporal random walk starting from lat_center, lon_center

    Parameters
    ----------
    data : xarray.DataArray
        The xarray containing the data.
    lat_center : int
        Center latitude index.
    lon_center : int
        Center longitude index.
    initial_time_layer : int
        Initial time layer index.
    steps : int
        Number of steps in a self-avoiding random walk used to sample 
        the neighborhood in space and time of a given observation.
    p_back: float
        Probability of stepping back in time during the random walk 
        (default is 0.2).
    connectivity : int
        Number of  max spatial neighbors to consider (default is 10).
    avoid_nan: bool
        Wheter to avoid nans while sampling the neighborhood or not.
        (Default is True.)
    d : int, optional
        Embedding dimension used in the sampled series to estimate the
        permutation patterns (default is 3).
    tau : int, optional
        Time delay used in the sampled series to estimate the
        permutation patterns (default is 10).

    Returns
    -------
    tuple
        array containing timeseries sampled with a spatio-temporal random walk.
    """
    da                          = data
    time_dim, lat_dim, lon_dim  = "time", "latitude", "longitude"
    T, H, W                     = da.sizes[time_dim], da.sizes[lat_dim], da.sizes[lon_dim]
    rng                         = np.random.default_rng()


    # Track visited spacetime cells
    def in_bounds(i, j): return 0 <= i < H and 0 <= j < W

    # Helper to fetch a single value (small .isel keeps compute granular on Dask)
    def get_val(tt, ii, jj):
        #return da.isel({time_dim: tt, lat_dim: ii, lon_dim: jj}).values()
        val = da.isel({time_dim: tt, lat_dim: ii, lon_dim: jj}).to_array()
        # Convert 0-d numpy array → scalar
        return float(val) if np.ndim(val) == 0 else val.item()
    
    values = []
    for _ in range(ensemble):
        visited = set()
        visited.add((initial_time_layer, lat_center, lon_center))
        tempvalues = [get_val(initial_time_layer, lat_center, lon_center)]
        t, i, j = initial_time_layer, lat_center, lon_center 
        for step_idx in range(steps-1):
            if step_idx == 0:
                tt = t+1
            else:
                tt += 1
            #if at this time all the dataset in nan insert a nan value and continue
            if  bool(da.CHL.isel({time_dim: tt}).isnull().all()):
            #if  bool(da.RRS412.isel({time_dim: tt}).isnull().all()):
                tempvalues.append(np.nan)
                visited.add((tt, i, j))
#               print(f"all nan at time {tt}, inserting nan value")
                continue
            # Reset rule: if step is multiple of tau*d go back to the initial coordinates
            reset=True
            if reset:
                if step_idx > 0 and step_idx % (tau*d) == 0:
                    i, j = lat_center, lon_center
                    v = get_val(tt, i, j)
                    if (not avoid_nan) or np.isfinite(v):
                        tempvalues.append(v)
                        visited.add((tt, i, j))
                        continue  # skip normal random choice this step
                    else:
                        pass # continue with normal random choice if the center is nan
            space_cand = []
            radius = 1
            max_radius = min(connectivity,max(H, W)) # max search radius to avoid infinite loops
            while not space_cand and radius <= max_radius:  # expand until we find something or exhaust grid
            # Build all offsets at this radius
                offsets = [(di, dj) for di in range(-radius, radius+1) 
                        for dj in range(-radius, radius+1)
                           ]
                for di, dj in rng.permutation(offsets):
                    ii, jj = i + di, j + dj
                    if in_bounds(ii, jj):  # and (t, ii, jj) not in visited:
                        v = get_val(tt, ii, jj)
                        if (not avoid_nan) or np.isfinite(v):
                            space_cand.append((tt, ii, jj, v))

                radius += 1


            #if cannot find valid space_cand add nan: # may be change to return nan?
            if not space_cand:
                print(f"cannot find valid candidates, adding a nan at step {step_idx}")
                space_cand = [(tt, i, j, np.nan)]
                #break  # stuck
            k = rng.integers(len(space_cand))
            t, i, j, v = space_cand[k]

            # add visited point
            visited.add((t, i, j))
            tempvalues.append(v)
        
        values.append(tempvalues)
    
    return values


    

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
    local_arr = np.empty((n_local, 3), dtype=float) #size of indicators 

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
