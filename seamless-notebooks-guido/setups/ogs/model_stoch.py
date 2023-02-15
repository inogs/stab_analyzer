import os,sys
import numpy as np
import scipy.integrate
import pyfabm
import netCDF4 as nc
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
import math
import subprocess
#import statsmodels.api as sm
from scipy.stats import variation
import matplotlib.pyplot as plt
import netCDF4 as nc
from scipy import signal
from scipy.stats import pearsonr

import copy
def signalPower(x):
    return np.mean(x**2)
def SNRsystem(inputSig, outputSig):
    noise = outputSig-inputSig

    powS = signalPower(outputSig)
    powN = signalPower(noise)
    return 10*np.log10((powS-powN)/powN)


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

ensamble_counter=int(np.loadtxt('ensamble_counter.txt'))

command='/g100_scratch/userexternal/gocchipi/BFM_STO/' + str(ensamble_counter).zfill(6) + '/'
#output of the model for all samples
#ntrials = 100
ntrials = 10  
#D = 0.001  # Standard deviation.

t_eval = np.linspace(0, 3650.*86400, 50000) 
y = np.zeros((ntrials,len(t_eval),54))
fig, ax = plt.subplots(3,3)
ax = ax.ravel()
fig1, ax1 = plt.subplots(3,3)
ax1 = ax1.ravel()
k=0

##########sinusoidal temperature fluctuation
# parameters
MeanTemp = 15           # Average temperature in the country
DailyAmpl = 5           # Amplitude of the daily cycle
YearlyAmpl = 5           # Amplitude of the yearly cycle
#oiseStd = 0.1          # Standard deviation of normally distributed error term

# Total seconds in year
TotalHours = 24*365*60*60 #year period
tau        = 24*60*60     #day period

# Generate the frequency components of the data
DailyCycle = -DailyAmpl*np.cos( (2*np.pi)*t_eval/tau )
YearlyCycle = -YearlyAmpl*np.cos( (2*np.pi)*t_eval/TotalHours )
#Noise = np.random.normal(0, NoiseStd, TotalHours)

# Final series
T = MeanTemp + DailyCycle + YearlyCycle #+ Noise




##########sinusoidal light fluctuation
# parameters
MeanTemp =  0           # Average temperature in the country
DailyAmpl = 100          # Amplitude of the daily cycle
YearlyAmpl = 10         # Amplitude of the yearly cycle
#NoiseStd = 0.1          # Standard deviation of normally distributed error term

# Total seconds in year
TotalHours = 24*365*60*60 #year period
tau        = 24*60*60     #day period

# Generate the frequency components of the data
DailyCycle = -DailyAmpl*np.cos( (2*np.pi)*t_eval/tau )
YearlyCycle = -YearlyAmpl*np.cos( (2*np.pi)*t_eval/TotalHours )
#Noise = np.random.normal(0, NoiseStd, TotalHours)

# Final series
L = np.where(MeanTemp + DailyCycle + YearlyCycle > 0, MeanTemp + DailyCycle + YearlyCycle, 0.)#+ Noise


#subprocess.run(["mkdir","-p",command])
#subprocess.run(["cp","fabm.yaml",command])

for s in range(ntrials):
    r = s % nranks
    if r == rank :
      # Create model (loads fabm.yaml)
      model = pyfabm.Model('fabm.yaml')
      
      # Configure the environment
      # Note: the set of environmental dependencies depends on the loaded biogeochemical model.
      model.dependencies['cell_thickness'].value = 1.
      model.dependencies['temperature'].value = T[0]#15.
      model.dependencies['practical_salinity'].value = 30.
      model.dependencies['density'].value = 1000.
      model.dependencies['depth'].value = 1.
      model.dependencies['pressure'].value = 1.
      model.dependencies['isBen'].value = 1.
      model.dependencies['longitude'].value = 0.
      model.dependencies['latitude'].value = 0.
      model.dependencies['surface_downwelling_shortwave_flux'].value = L[0]#50.
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
#     varnames = ["Z3/c"]
      varnames = ["B1/c","P1/c","P2/c","P3/c","P4/c","Z3/c","Z4/c","Z5/c","Z6/c"]

      #noise std deviation

      for i,v in enumerate(model.parameters):
         if 'D_noise' in v.name:
              D = v.getValue()
      print('noise std ', D)
      
      # Time-integrate over 1000 days (note: FABM's internal time unit is seconds!)
      dt = (t_eval[-1]-t_eval[0])/len(t_eval)
      sqrtdt = np.sqrt(dt)
      y[s,0,:]=model.state[:]
      
      #print(model.state) 
      for i in range(len(t_eval)):
          model.dependencies['cell_thickness'].value = 1.
          model.dependencies['temperature'].value = T[i]
          model.dependencies['practical_salinity'].value = 30.
          model.dependencies['density'].value = 1000.
          model.dependencies['depth'].value = 1.
          model.dependencies['pressure'].value = 1.
          model.dependencies['isBen'].value = 1.
          model.dependencies['longitude'].value = 0.
          model.dependencies['latitude'].value = 0.
          model.dependencies['surface_downwelling_shortwave_flux'].value = L[i]
          model.dependencies['surface_air_pressure'].value = 1.
          model.dependencies['wind_speed'].value = 5.
          model.dependencies['mole_fraction_of_carbon_dioxide_in_air'].value = 390.
          model.dependencies['number_of_days_since_start_of_the_year'].value = 1.
          dy = model.getRates()
          for j in range(len(model.state)):
              if i!=0:
                  if model.state_variables[j].name in varnames:
                      if i<15000:    #stop noise
                          y[s,i,j]=y[s,i-1,j]+dy[j]*dt  + D * sqrtdt * np.random.randn() * y[s,i-1,j]
                      else:
                          y[s,i,j]=y[s,i-1,j]+dy[j]*dt
                  else:
                      y[s,i,j]=y[s,i-1,j]+dy[j]*dt
          model.state[:]=y[s,i,:]
if rank == 0:
    y_global = np.copy(y)
    for s in range(ntrials):
        r = s % nranks
        if r !=0 :
          rea = comm.recv(source=r, tag=11)
          y_global[int(rea['idx']),:,:] = rea['data']
        comm.Barrier()
else:
    for s in range(ntrials):
        r = s % nranks
        rea = {'idx':s, 'data': y[s,:,:]}
        if r == rank :
            comm.send(rea, dest=0, tag=11)
        comm.Barrier()
comm.Barrier()
print('End of the loop',flush=True)

if rank == 0 :
  t_eval = t_eval / 86400
  y = np.copy(y_global)  
  # We create bins for the histograms.
  bins = np.linspace(0., 50., 100)
  
  #get deterministic trajectory
  ncname = 'result_ode.nc'
  f = nc.Dataset(ncname)
  varnames = ["B1/c","P1/c","P2/c","P3/c","P4/c","Z3/c","Z4/c","Z5/c","Z6/c"]
  varnames_ = []
  for v in varnames:
   varnames_.append(v.replace("/","_"))
  z = np.zeros((len(varnames),len(t_eval)))
  for it,tname in enumerate(varnames_):  
            z[it] = f.variables[tname][:,0,0,0]


  
  # Save results
  
  
# fileoutput = 'result.nc'
  fileoutput = 'result_4.nc'
  f = nc.Dataset(fileoutput, mode='w')
  
  lat_dim = f.createDimension('lat', 1)
  lon_dim = f.createDimension('lon', 1)
  dep_dim = f.createDimension('z', 1)
  time_dim = f.createDimension('time', t_eval.shape[0])
  #n=np.arange(0,ntrials,1)
  #sample_dim = f.createDimension('n', n.shape[0])
  
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
  f.variables['time'][:]=t_eval
  
  depth = f.createVariable('z', np.float32, ('z',))
  depth.units = 'meters'
  depth.long_name = 'depth'
  f.variables['z'][:]=1
  
  temp = f.createVariable('temp', np.float64, ('time',))
  temp.units = 'C'
  temp.long_name = 'temperature in celsius'
  f.variables['temp'][:]=T
  #sample = f.createVariable('n', np.float64, ('z',))
  #sample.units = ''
  #sample.long_name = 'sample'
  #f.variables['n'][:]=n
  
  for v,variable in enumerate(model.state_variables):
     ncvar = variable.name.replace("/","_")
     var = f.createVariable(ncvar, np.float64, ('time', 'z','lat','lon'))
     var.units = variable.units
     var.long_name = variable.long_name
     f.variables[ncvar][:]=np.mean(y[:,:,v], axis=0)
  #
  # sensitivity targets
  laststeps = int(len(t_eval)/2.)
  # cv
  y = np.mean(y,axis=0)
  kkk=0
  for v,variable in enumerate(model.state_variables):
     ncvar = variable.name.replace("/","_")
     if ncvar in varnames_:
         cv = (variation(y[-laststeps:,v])-variation(z[kkk,-laststeps:]))/D
  #      cv = pearsonr(y[-laststeps:,v], z[kkk,-laststeps:])[0]
  #      mean = np.mean(y[-laststeps:,v])/np.mean(z[kkk,-laststeps:])
         mean = np.mean(z[kkk,-laststeps:-100]*y[-(laststeps-100):,v])
         SNR = SNRsystem(y[-laststeps:,v], z[kkk,-laststeps:])
         if y[:,v].min() < 0.0001:
             ext=1
         else:
             ext=0
  #fourier sectrum
  #  if ncvar == 'Z6_c':
  #      plt.clf()
  #      PF = fft(y[-laststeps:,v])
  #      deltaT=t_eval[1]-t_eval[0]
  #      freq = (1/deltaT) * np.linspace(0,laststeps/2,int(laststeps/2)) / laststeps
  #      PF2 = np.abs(PF / laststeps) #full spectrum  /Npoints to get the true amplitudes
  #      PF1 = PF2[:laststeps//2]              #half spectrum
  #      PF1[1:] = 2*PF1[1:]
  #      plt.plot(freq,PF1)
  #      plt.show()
         
         ncfstate = variable.name.replace("/","_") + "_cv" 
         P1_fstate = f.createVariable(ncfstate, np.float32, ('z','lat','lon'))
         P1_fstate.units = ''
         P1_fstate.long_name = variable.name.replace("/","_") + "_cv"
         f.variables[ncfstate][:]=cv
         print(ncfstate,cv)
  
         ncfstate = variable.name.replace("/","_") + "_mean"
         P1_fstate = f.createVariable(ncfstate, np.float32, ('z','lat','lon'))
         P1_fstate.units = ''
         P1_fstate.long_name = variable.name.replace("/","_") + "_mean"
         f.variables[ncfstate][:]=mean
         print(ncfstate,mean)
  
         ncfstate = variable.name.replace("/","_") + "_snr"
         P1_fstate = f.createVariable(ncfstate, np.float32, ('z','lat','lon'))
         P1_fstate.units = ''
         P1_fstate.long_name = variable.name.replace("/","_") + "_snr"
         f.variables[ncfstate][:]=SNR
         print(ncfstate,SNR)
         print('*****')
  
         ncfstate = variable.name.replace("/","_") + "_ext"
         P1_fstate = f.createVariable(ncfstate, np.float32, ('z','lat','lon'))
         P1_fstate.units = ''
         P1_fstate.long_name = variable.name.replace("/","_") + "_ext"
         f.variables[ncfstate][:]=ext
  
  f.close()
#subprocess.run(["cp","result.nc",command])
