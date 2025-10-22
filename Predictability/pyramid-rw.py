from mpi4py import MPI
import numpy as np
import xarray as xr
import pandas as pd
import lyapunovV
#import nolds as NOLDS

from os import cpu_count
from matplotlib import pyplot as plt
import ordpy as od

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


def SAMPLED_HC_LAT_LON(data, lat_center, lon_center, initial_time_layer, steps=10, ensemble=10, p_back=0.2, connectivity=10, 
                       avoid_nan=True, d=3, tau=10):
    """
    Calculate permutation complexity and permutation entropy for a given xarray dataset using random-walk sampling
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
        permutation entropy, permutation complexity, valid sample size, lyapunov exponent, nolds lyapunov

    """
    values = SAMPLED_LAT_LON(data, lat_center, lon_center, initial_time_layer, steps, ensemble, p_back, connectivity, 
                       avoid_nan, d, tau)
    
    #check if values is empty
    if len(values) == 0:
        return np.nan, np.nan, 0, np.nan, np.nan    
    #compute PE, C, lyapunov exponent and nolds on the axis=1 of values and than average
    PE_list= np.empty(len(values), dtype=float)
    C_list  = np.empty(len(values), dtype=float)
    lyap_list = np.empty(len(values), dtype=float)
    nolds_list = np.empty(len(values), dtype=float)
    fig,ax = plt.subplots()
    for i in range(len(values)):
        ts = np.array(values[i][:])

        ### WATCH OUT
        
        #remove nans to compute indicators, in the dataset there is some timesetp where all is nan
        # should be removed in final implementation
        valid_ts = ts[~np.isnan(ts)]

        ####
        #when removed use
    #    valid_ts = ts
        valid_sample_size = len(valid_ts)
        
        #plot values to check the random walk results, can be removed later
        valid_ts = valid_ts[-2000:]
        ax.plot(valid_ts,alpha=0.3)
        fig.savefig('prova3.png')

        #check if valid sample size is enough
        if valid_sample_size < steps*0.3:
            lyap_list[i] = np.nan
            nolds_list[i] = np.nan
            PE_list[i] = np.nan
            C_list[i] = np.nan
        else:
            #compute PE, C, lyapunov exponent and nolds
            try:
                PE_list[i], C_list[i] = od.complexity_entropy(valid_ts, dx=d, taux=tau)
            except:
                PE_list[i], C_list[i] = np.nan, np.nan
            lp = lyapunovV.LYAP(valid_ts)
            try:
                lyap_list[i] = lp.lyap_e(dt=1,ndim=3,evolve=20,tau=tau)
            except:
                lyap_list[i] = np.nan
            try:
                #remove nolds calculation for the time being
                nolds_list[i] = np.nan
                #nolds_list[i] = NOLDS.lyap_r(valid_ts, emb_dim=3, tau=tau)
            except:
                nolds_list[i] = np.nan
    PE = np.nanmean(PE_list)
    C  = np.nanmean(C_list)
    lyap = np.nanmean(lyap_list)
    nolds = np.nanmean(nolds_list)

    return PE, C, len(values), lyap, nolds
    

#used to run over coordinates
def mpi_pyramidal_complexity_entropy_map(
    data, initial_time_layer, delta_time, n_blocks_time,
    dx=2, dy=2, dz=2, N_ens=10, rank=0, size=1
):
    """
    MPI-parallel version (mpi4py) of the original Dask implementation.
    - On rank 0, returns an xr.Dataset with vars: 'entropy', 'complexity', 'valid_sample_size'.
      On non-root ranks, returns None.
    Calculate permutation complexity and permutation entropy map for a given xarray dataset using random-walk sampling

    Parameters
    ----------
    data : xarray.Dataset
        The xarray dataset containing the data.
    initial_time_layer : int
        Initial time layer index.
    delta_time : int
        Time steps for blocks.
    n_blocks_time : int
        Number of time blocks to consider.
    dx : int, optional
        Embedding dimension in the x (longitude) direction (default is 2).
    dy : int, optional
        Embedding dimension in the y (latitude) direction (default is 2).
    dz : int, optional
        Embedding dimension in the z (time) direction (default is 2).
    N_ens : int
        Number of ensemble members to run a random walk for each point (default is 10).

    Returns
    -------
    xr.Dataset with coords ('latitude','longitude') and data_vars
        'entropy', 'complexity', 'valid_sample_size'
    """

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

    # Preallocate local results (columns: idx, entropy, complexity, valid)
    n_local = len(my_indices)
    local_arr = np.empty((n_local, 6), dtype=float)

    for i, idx in enumerate(my_indices):
        lat = idx // n_lon
        lon = idx % n_lon
        ent, comp, valid, lyap, nolds = SAMPLED_HC_LAT_LON(
            data, lat, lon,
            initial_time_layer,
            tau=delta_time,
            steps=n_blocks_time, ensemble=N_ens, 
            p_back=0.5, connectivity=8, 
            avoid_nan=True, d=dx)
        local_arr[i, 0] = idx
        local_arr[i, 1] = ent
        local_arr[i, 2] = comp
        local_arr[i, 3] = valid
        local_arr[i, 4] = lyap
        local_arr[i, 5] = nolds

    # Gather arrays at root
    gathered = comm.gather(local_arr, root=0)

    if rank == 0:
        # Stack all results
        all_results = np.vstack(gathered)

        # Allocate flat arrays
        ent_flat   = np.empty(total, dtype=float)
        comp_flat  = np.empty(total, dtype=float)
        valid_flat = np.empty(total, dtype=float)
        lyap_flat  = np.empty(total, dtype=float)
        nolds_flat = np.empty(total, dtype=float)

        # Fill arrays using indices
        idxs = all_results[:, 0].astype(int)
        ent_flat[idxs]   = all_results[:, 1]
        comp_flat[idxs]  = all_results[:, 2]
        valid_flat[idxs] = all_results[:, 3]
        lyap_flat[idxs] = all_results[:, 4]
        nolds_flat[idxs] = all_results[:, 5]

        # Reshape into (lat, lon)
        ent   = ent_flat.reshape((n_lat, n_lon))
        comp  = comp_flat.reshape((n_lat, n_lon))
        valid = valid_flat.reshape((n_lat, n_lon))
        lyap  = lyap_flat.reshape((n_lat, n_lon))
        nolds = nolds_flat.reshape((n_lat, n_lon))

        return xr.Dataset(
            {
                'entropy':           (('latitude','longitude'), ent),
                'complexity':        (('latitude','longitude'), comp),
                'valid_sample_size': (('latitude','longitude'), valid),
                'lyapunov':          (('latitude','longitude'), lyap),
                'nolds':             (('latitude','longitude'), nolds),
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

# Load data regional test dataset

data = xr.open_dataset('/g100_scratch/userexternal/gocchipi/GLOBAL/tmp_daily/regional_dataset_CHL_10_20_-30_-15.nc', decode_times=False)
#data = xr.open_dataset('/g100_scratch/userexternal/gocchipi/DATA_for_PE/regional_dataset_CHL_10_20_-30_-15.nc', decode_times=False)
#data = xr.open_dataset('/g100_scratch/userexternal/gocchipi/DATA_for_PE/regional_dataset_RRS412_10_20_-30_-15.nc', decode_times=False)

# select a smaller region for rapid testing
data = data.sel(latitude=slice(17, 18), longitude=slice(-18, -17))

units, reference_date = data.time.attrs['units'].split('since')

#convert time to datetime64
data['time'] = pd.date_range(start=reference_date, periods=data.sizes['time'], freq='d')

#clean data removing negative and too high values
data['CHL'] = data['CHL'].where((data['CHL'] > 0) & (data['CHL'] < 100))

# choose parameters for testing
array_delta_time = [1,10,30,180,360]  #time delays to test (tau)
array_N_ens = [5,10,30]               #number of ensemble members to test

fig,axs = plt.subplots(4, 15, figsize=(20, 50))

for ii, delta_time in enumerate(array_delta_time):
    for jj, N_ens in enumerate(array_N_ens):
        print(f"Run delta_time: {delta_time}, N_ens: {N_ens}",flush=True)
        fraction = 1.  # use 100% of data, can choose smaller than one to reduce the computation time
        initial_time_layer = 0
        n_blocks_time = int((data.sizes['time'] - initial_time_layer) * fraction)

        # run parallel computation
        res_data = mpi_pyramidal_complexity_entropy_map(data, initial_time_layer=initial_time_layer,
                                 delta_time=delta_time, n_blocks_time=n_blocks_time, dx=6,N_ens=N_ens)
        # plot and save results
        if rank ==0:
            #plot results
            res_data.lyapunov.plot(ax=axs[0,int(ii*len(array_N_ens)+jj)], cmap='viridis', vmin=np.nanmin(res_data.lyapunov), vmax=np.nanmax(res_data.lyapunov)+1.e-2)
            res_data.nolds.plot(ax=axs[1,int(ii*len(array_N_ens)+jj)], cmap='viridis', vmin=np.nanmin(res_data.nolds), vmax=np.nanmax(res_data.nolds)+1.e-2)
            res_data.entropy.plot(ax=axs[2,int(ii*len(array_N_ens)+jj)], cmap='viridis', vmin=np.nanmin(res_data.entropy), vmax=np.nanmax(res_data.entropy)+1.e-2)
            res_data.complexity.plot(ax=axs[3,int(ii*len(array_N_ens)+jj)], cmap='viridis', vmin=np.nanmin(res_data.complexity), vmax=np.nanmax(res_data.complexity)+1.e-2)

            axs[0,int(ii*len(array_N_ens)+jj)].set_title(f"Lyap tau: {delta_time}, N_ens: {N_ens}")
            axs[1,int(ii*len(array_N_ens)+jj)].set_title(f"Nolds tau: {delta_time}, N_ens: {N_ens}")
            axs[2,int(ii*len(array_N_ens)+jj)].set_title(f"PE tau: {delta_time}, N_ens: {N_ens}")
            axs[3,int(ii*len(array_N_ens)+jj)].set_title(f"C tau: {delta_time}, N_ens: {N_ens}")
            #intermediate save to runtime checks
            fig.savefig('RW_indicators.png', dpi=300)
            
            #save results
            res_data.to_netcdf(f'RW_indicators_{delta_time}_{N_ens}.nc', mode='w', format='NETCDF4')

fig.tight_layout()
fig.savefig('RW_indicators.png', dpi=300)

