#! /bin/bash

#SBATCH -N 2
#SBATCH --ntasks-per-node=48
#SBATCH --time=04:00:00
#SBATCH --mem=300gb
#SBATCH --job-name=pyramid-PE
#SBATCH --account OGS23_PRACE_IT
#SBATCH --partition g100_usr_prod
##SBATCH --qos=g100_qos_dbg

export MPI4PY_FUTURES_MAX_WORKERS=48
eval "$(/g100/home/userexternal/gocchipi/my_libraries/conda-2023/bin/conda shell.bash hook)"
conda activate conda-sto
module load intel

date
mpirun -np 48 python pyramid-rw.py
date
