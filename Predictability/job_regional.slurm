#! /bin/bash

#SBATCH -N 3
#SBATCH --ntasks-per-node=48
#SBATCH --time=12:00:00
#SBATCH --mem=300gb
#SBATCH --job-name=lat_Trends
#SBATCH --account OGS23_PRACE_IT
#SBATCH --partition g100_usr_prod
##SBATCH --qos=g100_qos_dbg

eval "$(/g100/home/userexternal/gocchipi/my_libraries/conda-2023/bin/conda shell.bash hook)"
#conda activate copernicus
conda activate conda-sto

module load intel

date

#python create_regional_dataset.py
mpirun -np 144 python regional_study.py

date
