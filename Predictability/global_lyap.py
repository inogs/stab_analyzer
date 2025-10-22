from mpi4py import MPI
import xarray as xr
import numpy as np
import lyapunovV
import glob
import ordpy as od
import os

import warnings
warnings.filterwarnings('ignore')

#use argparse to get two integer numbers the first passing the job number the second the total number of jobs
import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    read job id and total number of jobs from command line arguments.
    ''',
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(   '--jobid', '-id',
                                type = int,
                                required = True,
                                default = '',
                                help = 'number of the job to run, starting from 0')
    parser.add_argument(   '--totaljobs', '-tot',
                                type = int,
                                required = True,
                                default = '',
                                help = 'total number of jobs to run')
    return parser.parse_args()

args = argument()
id_job = args.jobid
total_jobs = args.totaljobs



# MPI setup
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# choose variable to process
var = 'CHL'
#var = 'RRS412'

# Get list of data files, check the path
if var == 'CHL':
    data_files = glob.glob('/g100_scratch/userexternal/gocchipi/GLOBAL/tmp_daily/cmems_????-*.nc')
elif var == 'RRS412':
    data_files = glob.glob('/g100_scratch/userexternal/gocchipi/DATA_for_PE/ALL/cmmems_RRS412_????-*.nc')
data_files.sort()

# Get coordinates
ds = xr.open_dataset(data_files[0], decode_times=False)
latitudes = ds.latitude.values
longitudes = ds.longitude.values
coords = np.array(np.meshgrid(latitudes, longitudes)).T.reshape(-1, 2)
ds = None

#Get a fraction of coordinates based on the job number and total number of jobs
coords = coords[id_job::total_jobs]


# Downsample to speed up
#oords = coords[::10000]
chunk_size = len(coords) // size
start = rank * chunk_size
end = None if rank == size - 1 else (rank + 1) * chunk_size
coords_chunk = coords[start:end]

# Save each rank output to a separate CSV file
output_file = f'CSV/lyapunov_{var}_rank{rank}.csv'
merged_file = 'CSV/merged_lyapunov.csv'

for ii, coord in enumerate(coords_chunk):
    try:
        # Check if the coordinate already exists in the merged file, if true skip the computation
        if os.path.exists(merged_file):
            with open(merged_file, 'r') as f:
                if any(f"{coord[0]},{coord[1]}" in line for line in f):
                    continue
        if os.path.exists(output_file):
            with open(output_file, 'r') as f:
                if any(f"{coord[0]},{coord[1]}" in line for line in f):
                    continue
    except:
        pass

    y = np.zeros(len(data_files))
    for jj, file in enumerate(data_files):
        ds = xr.open_dataset(file, decode_times=False)
        try:
            if var == 'CHL':
                # if ds.sel(latitude=coord[0], longitude=coord[1], method='nearest').CHL.values is not a number compute the mean of the array
                if np.shape(ds.sel(latitude=coord[0], longitude=coord[1], method='nearest').CHL.values) == ():
                    y[jj] = ds.sel(latitude=coord[0], longitude=coord[1], method='nearest').CHL.values
                else:
                    y[jj] = np.nanmean(ds.sel(latitude=coord[0], longitude=coord[1], method='nearest').CHL.values)
            elif var == 'RRS412':
                if np.shape(ds.sel(latitude=coord[0], longitude=coord[1], method='nearest').RRS412.values) == ():
                    y[jj] = ds.sel(latitude=coord[0], longitude=coord[1], method='nearest').RRS412.values
                else:
                    y[jj] = np.nanmean(ds.sel(latitude=coord[0], longitude=coord[1], method='nearest').RRS412.values)
        except:
            y[jj] = np.nan

    # remove invalid values
    y[(y < 0) | (y > 100)] = np.nan

    ### WATCH OUT ###
    # here remove nan values otherwise lyapunov do not works, but need to overcome the problem
    
    y = y[~np.isnan(y)]

    ######

    # skip if not enough data
    if len(y) < 10:
        continue
    print(f"Processing coordinates {coord} with {len(y)} valid data points",flush=True)
    try:
        lyap = lyapunovV.LYAP(y)
        lyap_val = lyap.lyap_e(dt=1)
    except:
        lyap_val = np.nan

    try:
        PE_val, C_val = od.complexity_entropy(y, dx=6)
    except:
        PE_val, C_val = np.nan, np.nan

    with open(output_file, 'a') as f:
        f.write(f"{coord[0]},{coord[1]},{lyap_val},{PE_val},{C_val}\n")

# Merge all CSV files from all ranks
if rank == 0:
    for r in range(size):
        rank_file = f'CSV/lyapunov_{var}_rank{r}.csv'
        if os.path.exists(rank_file):
            with open(rank_file, 'r') as rf:
                with open(merged_file, 'a') as mf:
                    for line in rf:
                        mf.write(line)

#save the results in a netcdf file with the same strucutre of ds but without time dimension
ds_out = xr.Dataset(
    {
        'lyapunov': (['latitude', 'longitude'], lyap_coords.reshape(ds.latitude.size, ds.longitude.size)),
        'PE': (['latitude', 'longitude'], PE.reshape(ds.latitude.size, ds.longitude.size)),
        'C': (['latitude', 'longitude'], C.reshape(ds.latitude.size, ds.longitude.size)),
    },
    coords={
        'latitude': ds.latitude,
        'longitude': ds.longitude,
    }
)
ds_out.to_netcdf(f'lyapunov_{var}.nc')
