# Scripts to compute predictability of ocean color 

## Timeseries analysis of each gridpoint

### Regional Test 
1. create_regional_dataset.py : extract a subset of global  data
2. regional_study.py : compute the metrics in the regional subset for different parameters (tau) values, better to use a parallel job (e.g. 2 nodes), needs folder REGIONAL to save outputs
3. cat REGIONAL/lyapunov_*.csv >> REGIONAL/merged.csv : merge multiple ranks outputs of previous step
4. plot_lyapunov_test.py: plot of metrics in merged.csv
5. job_regional.slurm : example of job to run create_regional_dataset.py and regional_study.py

### Global computation 
1. global_lyap.py : parallel script to compute metrics lyap, PE, C over the globe, needs folder CSV to save outputs
2. job_ensemble.slurm : example of job to run global_lyap.py
3. launcher.sh : example of launcher of array of jobs
4. cat CSV/lyapunov_*.csv >> CSV/merged_lyapunov.csv : merge multiple ranks outputs of previous step
5. plot_global_lyap.py: plot the result of global computation

## Spatio-temporal analysis with random-walk 
1. pyramid-rw.py : contains function to do random-walk analysis of spatio-temporal dataset, plus example of running over a regional subset for testing
2. job_rw.slurm : job to launch pyramid-rw.py in parallel
