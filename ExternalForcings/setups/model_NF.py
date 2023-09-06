import os,sys
import numpy as np
import scipy.integrate
import pyfabm
import netCDF4 as nc
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
import math
import subprocess
import matplotlib.pyplot as plt

ensamble_counter=int(np.loadtxt('ensamble_counter.txt'))

command='/g100_work/OGS21_PRACE_P/gocchipinti/PeriodicForcing/seamless-notebooks-guido/parsac/OUTPUTS/PF/' + str(ensamble_counter).zfill(6) + '/'


subprocess.run(["mkdir","-p",command])
subprocess.run(["cp","fabm.yaml",command])


try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False

print("Init ...", flush=True)
print(isParallel, flush=True)

#output of the model for all samples
#ntrials = 1000
ntrials = 1   
#t_eval = np.linspace(0, 3650.*86400, 50000) 
t_eval = np.linspace(0, 7300.*86400, 100000) 
y = np.zeros((ntrials,len(t_eval),54))
k=0
l=[]


#sinusoidal temperature fluctuation
# parameters
MeanTemp = 15           # Average temperature in the country
DailyAmpl = 0           # Amplitude of the daily cycle
YearlyAmpl =  5         # Amplitude of the yearly cycle
#oiseStd = 0.1          # Standard deviation of normally distributed error term

# Total seconds in year
TotalHours = 24*365*60*60 #year period
tau        = 24*60*60     #day period

# Generate the frequency components of the data
DailyCycle = -DailyAmpl*np.cos( (2*np.pi)*t_eval/tau )
YearlyCycle = -YearlyAmpl*np.cos( (2*np.pi)*t_eval/TotalHours )
#Noise = np.random.normal(0, NoiseStd, TotalHours)

# Final series
T = MeanTemp + YearlyCycle #+ Noise




varnames = ["B1/c","P1/c","P2/c","P3/c","P4/c","Z5/c","Z6/c","Z3/c","Z4/c"]

for s in range(ntrials):
    r= s % nranks
    if r == rank :
      # Create model (loads fabm.yaml)
      model = pyfabm.Model('fabm.yaml')
      
      # Configure the environment
      # Note: the set of environmental dependencies depends on the loaded biogeochemical model.
      model.dependencies['cell_thickness'].value = 1.
      model.dependencies['temperature'].value = T[0]#15.
#     model.dependencies['temperature'].value = 15.
      model.dependencies['practical_salinity'].value = 30.
      model.dependencies['density'].value = 1000.
      model.dependencies['depth'].value = 1.
      model.dependencies['pressure'].value = 1.
      model.dependencies['isBen'].value = 1.
      model.dependencies['longitude'].value = 0.
      model.dependencies['latitude'].value = 0.
      model.dependencies['surface_downwelling_shortwave_flux'].value = 50.
#     model.dependencies['surface_downwelling_shortwave_flux'].value = L[0]#50.
      model.dependencies['surface_air_pressure'].value = 1.
      model.dependencies['wind_speed'].value = 5.
      model.dependencies['mole_fraction_of_carbon_dioxide_in_air'].value = 390.
      model.dependencies['number_of_days_since_start_of_the_year'].value = 1.
      
      # Verify the model is ready to be used
      model.cell_thickness=1.
      
      assert model.checkReady(), 'One or more model dependencies have not been fulfilled.'
      
      # Time derivative
      #def dy(t0, y):
      #    model.state[:] = y
      #    return model.getRates()
      
      #variables with noise
      varnames = ["B1/c","P1/c","P2/c","P3/c","P4/c","Z5/c","Z6/c","Z3/c","Z4/c"]
      #noise
      D = 0  # Standard deviation.
      
      # Time-integrate over 1000 days (note: FABM's internal time unit is seconds!)
      dt = (t_eval[-1]-t_eval[0])/len(t_eval)
      sqrtdt = np.sqrt(dt)
      y[s,0,:]=model.state[:]
      
      #print(model.state) 
      for i in range(len(t_eval)):
          model.dependencies['temperature'].value = T[i]
#         model.dependencies['surface_downwelling_shortwave_flux'].value = L[i]
          dy = model.getRates()
          for j in range(len(model.state)):
              if i!=0:
                  y[s,i,j]=y[s,i-1,j]+dy[j]*dt
          model.state[:]=y[s,i,:]
      #if s!=ntrials-1:
      #    del model
if rank == 0:
    y_global = np.copy(y)
    val = np.zeros(y.shape)
    for s in range(ntrials):
        r = s % nranks
        if r !=0 :
          print(r)
          comm.Recv( [val[:], MPI.DOUBLE], source=r, tag=1 )
          y_global[s,:] = val[s,:]
else:
    val = np.zeros(y.shape)
    for s in range(ntrials):
        r = s % nranks
        if r == rank :
            comm.Send( [y[s,:], MPI.DOUBLE], dest=0, tag=1 )
print('End of the loop',flush=True)

if rank == 0 :
  y = np.copy(y_global)  
  
  
t = t_eval/86400

# Plot results
#pylab.plot(t, y)
#pylab.legend([variable.path for variable in model.state_variables])
#pylab.show()


Nt=t.shape[0]
deltaT=t[1]-t[0]
laststeps = int(Nt/2) #compute the indicators just for this last timesteps
freq = (1/deltaT) * np.linspace(0,laststeps/2,int(laststeps/2)) / laststeps

# Save results


fileoutput = 'result.nc'
f = nc.Dataset(fileoutput, mode='w')

lat_dim = f.createDimension('lat', 1)
lon_dim = f.createDimension('lon', 1)
dep_dim = f.createDimension('z', 1)
time_dim = f.createDimension('time', Nt)

lat = f.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
f.variables['lat'][:]=0

lon = f.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
f.variables['lon'][:]=0

time = f.createVariable('time', np.float64, ('time',))
time.units = 'days'
time.long_name = 'time'
f.variables['time'][:]=t

depth = f.createVariable('z', np.float32, ('z',))
depth.units = 'meters'
depth.long_name = 'depth'
f.variables['z'][:]=1

temp = f.createVariable('temp', np.float64, ('time',))
temp.units = 'C'
temp.long_name = 'temperature in celsius'
f.variables['temp'][:]=T


for v,variable in enumerate(model.state_variables):
  if variable.name in varnames:
   ncvar = variable.name.replace("/","_")
   var = f.createVariable(ncvar, np.float64, ('time', 'z','lat','lon'))
   var.units = variable.units
   var.long_name = variable.long_name
   f.variables[ncvar][:]=y[0,:,v]

## Sensibility targets
#
y = y[0,:,:]
#critical lenght for cycles
epsilon = 0.001
lyapunov = np.zeros(len(varnames))
ii=0
extinction = 0
for v,variable in enumerate(model.state_variables):
  if variable.name in varnames:
#coeff of variance
    p1_var = 0
    cv = np.var(y[-laststeps:,v])/np.mean(y[-laststeps:,v])
    if math.isinf(cv) :
        p1_var=float('nan')
#   elif variation(y[-laststeps:,v]) > epsilon :
#       p1_var = int(1)
    else :
        p1_var = cv

#Relative peaks in the fourier spac
    PF = fft(y[-laststeps:,v])

    PF2 = np.abs(PF / laststeps) #full spectrum  /Npoints to get the true amplitudes
    PF1 = PF2[:laststeps//2]              #half spectrum
    PF1[1:] = 2*PF1[1:]           #actual amplitude

    #find peaks
    peak_index0, dict_vals = find_peaks(np.concatenate(([min(PF1)],PF1,[min(PF1)])), distance=5)
    #find the biggest two
    peak_index0 = peak_index0-1
    heights0 = np.sort(PF1[peak_index0])[::-1]
    p1_peak = 0
    if len(heights0) >1 :
        if heights0[0] >0.001 :
            p1_peak = heights0[1]/heights0[0]
        else :
            p1_peak = 0
    else :
        p1_peak = 0
    #make it discrete
    if p1_peak > epsilon :
        p1_peak = int(1)
    else :
        p1_peak = int(0)

    p1_var = p1_var * p1_peak
    
    ncvar = variable.name.replace("/","_") + "_var"
    P1_var = f.createVariable(ncvar, np.float32, ('z','lat','lon'))
    P1_var.units = ' '
    P1_var.long_name = ncvar + 'iance'
    f.variables[ncvar][:]=p1_var

    #extinction
    
    if y[-laststeps:,v].max() < 1.e-4:
        extinction = 1



    #Lyapunov exponents
    if p1_var == 0:
        continue
    else:
              diffvar = y[-50000:,v]
              initial1 = 0
              initial2 = 5000
              diff = np.zeros(5000)
              for kk in range(int(len(diffvar)/5000-2)):
                      diff += diffvar[initial2:initial2+5000] - diffvar[initial1:initial1+5000]
                      initial1 += 5000
                      initial2 += 5000
              diff = diff/np.mean(diffvar)#/int(len(diffvar/5000)-3)
              if np.max(diff) < 1.e-3:
                  lyapunov[ii] = 0
              else:
                  import lyapunovV
                  lyap = lyapunovV.LYAP(y[-20000:,v])
                  try:
                    lyapunov[ii] = lyap.lyap_e(dt=7300./100000.)
                  except:
                    lyapunov[ii] = 0
              ii += 1


max_lyap = lyapunov.max()

ncvar = "max_lyapunov"
P1_var = f.createVariable(ncvar, np.float32, ('z','lat','lon'))
P1_var.units = ' '
P1_var.long_name = 'maximum_lyapunov_exponent'
f.variables[ncvar][:]=max_lyap

ncvar = "extinction"
P1_var = f.createVariable(ncvar, np.float32, ('z','lat','lon'))
P1_var.units = ' '
P1_var.long_name = 'extinction'
f.variables[ncvar][:]=extinction


f.close()
subprocess.run(["cp","result.nc",command])
