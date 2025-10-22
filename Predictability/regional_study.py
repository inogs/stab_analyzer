from mpi4py import MPI
import xarray as xr
import numpy as np
import lyapunovV
import glob
import ordpy as od
import os
import nolds

import warnings
warnings.filterwarnings('ignore')


# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Choose which variable to process
var = 'CHL'
#var = 'RRS412'

# Assign the correct data file path according to the variable
if var == 'CHL':
    data_file = '/g100_scratch/userexternal/gocchipi/DATA_for_PE/regional_dataset_CHL_10_20_-30_-15.nc'
elif var == 'RRS412':
    data_file = '/g100_scratch/userexternal/gocchipi/DATA_for_PE/regional_dataset_RRS412_10_20_-30_-15.nc'

# Get coordinates
ds = xr.open_dataset(data_file, decode_times=False)
# Restrict latitudes and longitudes if need to process a smaller region
ds = ds.sel(latitude=slice(16, 18), longitude=slice(-18, -15))
###
latitudes = ds.latitude.values
longitudes = ds.longitude.values
coords = np.array(np.meshgrid(latitudes, longitudes)).T.reshape(-1, 2)

#create empty arrays to store results
lyap_val = np.zeros(coords.shape[0])
PE_val= np.zeros(coords.shape[0])
C_val = np.zeros(coords.shape[0])



# Split coordinates among MPI ranks
chunk_size = len(coords) // size
start = rank * chunk_size
end = (rank + 1) * chunk_size if rank != size - 1 else len(coords)
coords_chunk = coords[start:end]

# Sace output file as a csv per rank to avoid write conflicts
output_file = f'REGIONAL/lyapunov_{var}_rank{rank}.csv'

#Compute the metrics
for ii, coord in enumerate(coords_chunk):
    try:
        # Skip if coordinate already in output file
        if os.path.exists(output_file):
            with open(output_file, 'r') as f:
                if any(f"{coord[0]},{coord[1]}" in line for line in f):
                    continue
    except:
        pass

    y = ds.sel(latitude=coord[0], longitude=coord[1], method='nearest')[var].values
   
    # remove invalid data
    y[(y < 0) | (y > 100)] = np.nan
#    y = y[~np.isnan(y)]
    
    # skip if not enough data points
    if len(y) < 1000:
        continue

    print(f"[Rank {rank}] Processing {coord} with {len(y)-np.sum(np.isnan(y))} valid points", flush=True)
    lp = lyapunovV.LYAP(y)
    #compute metrics for different tau values
    try:
        lyap1 = lp.lyap_e(dt=1,ndim=6,evolve=20,tau=1)
    except: 
        lyap1 = np.nan
    try:
        lyap10 = lp.lyap_e(dt=1,ndim=6,evolve=20,tau=10)
    except: 
        lyap10 = np.nan
    try:
        lyap30 = lp.lyap_e(dt=1,ndim=6,evolve=20,tau=30)
    except: 
        lyap30 = np.nan
    try:
        lyap180 = lp.lyap_e(dt=1,ndim=6,evolve=20,tau=180)
    except: 
        lyap180 = np.nan
    try:
        lyap360 = lp.lyap_e(dt=1,ndim=6,evolve=20,tau=360)
    except: 
        lyap360 = np.nan
    try:
        lyapnolds1 = nolds.lyap_r(y, emb_dim=6, tau=1)
    except:
        lyapnolds1 = np.nan
    try:
        lyapnolds10 = nolds.lyap_r(y, emb_dim=6, tau=10)
    except:
        lyapnolds10 = np.nan
    try:
        lyapnolds30 = nolds.lyap_r(y, emb_dim=6, tau=30)
    except:
        lyapnolds30 = np.nan
    try:
        lyapnolds180 = nolds.lyap_r(y, emb_dim=6, tau=180)
    except:
        lyapnolds180 = np.nan
    try:
        lyapnolds360 = nolds.lyap_r(y, emb_dim=6, tau=360)
    except:
        lyapnolds360 = np.nan

    try:
        pe1, c1 = od.complexity_entropy(y, dx=6, taux=1)
        pe10,c10 = od.complexity_entropy(y, dx=6, taux=10)
        pe30,c30 = od.complexity_entropy(y, dx=6, taux=30)
        pe180, c180 = od.complexity_entropy(y, dx=6, taux=180)
        pe360, c360 = od.complexity_entropy(y, dx=6, taux=360)
    except:
        pe1, c1 = np.nan, np.nan
        pe10, c10 = np.nan, np.nan
        pe30, c30 = np.nan, np.nan
        pe180, c180 = np.nan, np.nan
        pe360, c360 = np.nan, np.nan

    # Append results to the output file
    with open(output_file, 'a') as f:
        f.write(f"{coord[0]},{coord[1]},{lyap1},{lyap10},{lyap30},{lyapnolds1},{lyapnolds10},{lyapnolds30},{pe1},{c1},{pe10},{c10},{pe30},{c30},{lyap180},{lyap360},{lyapnolds180},{lyapnolds360},{pe180},{pe360},{c180},{c360}\n")
    

