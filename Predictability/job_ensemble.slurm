#! /bin/bash
#
#SBATCH --job-name=dem_noise
#SBATCH -N5                     
#SBATCH --ntasks-per-node=48    
#SBATCH --time=24:00:00
##SBATCH --mem=100gb
#SBATCH --account=OGS23_PRACE_IT
#SBATCH --partition=g100_usr_prod
#SBATCH --array=0-499            # Job array over trials (0 to 49)

echo $SLURM_JOB_ID
echo $SLURM_ARRAY_JOB_ID
echo $SLURM_ARRAY_TASK_ID
echo $SLURM_ARRAY_TASK_COUNT
echo $SLURM_ARRAY_TASK_MAX
echo $SLURM_ARRAY_TASK_MIN

cd $SLURM_SUBMIT_DIR

ulimit -s unlimited

date

cd $SLURM_SUBMIT_DIR
eval "$(/g100/home/userexternal/gocchipi/my_libraries/conda-2023/bin/conda shell.bash hook)"
conda activate conda-sto
module load intel
# Number of job array trials
# Need to choose based on the number of jobs in which subdivide the globe
# Must be consistent with the number of jobs in laucher.sh
export N_trials=500              


export i=${SLURM_ARRAY_TASK_ID} # Current trial index
export BASEDIR=$PWD

# Run the Python script with the current trial index and total number of trials 
mpirun -np 240 python global_lyap.py -id $i -tot $N_trials 


date