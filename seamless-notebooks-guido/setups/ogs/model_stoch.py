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


try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
    print(rank)
    print(nranks)
except:
    rank   = 0
    nranks = 1
    isParallel = False


#output of the model for all samples
ntrials = 1000
#ntrials = 100
#D = 0.001  # Standard deviation.

t_eval = np.linspace(0, 3650.*86400, 50000) 
y = np.zeros((ntrials,len(t_eval),54))
A = np.zeros((ntrials,len(t_eval),54))
fig, ax = plt.subplots(3,3)
ax = ax.ravel()
fig1, ax1 = plt.subplots(3,3)
ax1 = ax1.ravel()
k=0


for s in range(ntrials):
    r = s % nranks
    if r == rank :
      # Create model (loads fabm.yaml)
      model = pyfabm.Model('fabm.yaml')
      
      # Configure the environment
      # Note: the set of environmental dependencies depends on the loaded biogeochemical model.
      model.dependencies['cell_thickness'].value = 1.
      model.dependencies['temperature'].value = 15.
      model.dependencies['practical_salinity'].value = 30.
      model.dependencies['density'].value = 1000.
      model.dependencies['depth'].value = 1.
      model.dependencies['pressure'].value = 1.
      model.dependencies['isBen'].value = 1.
      model.dependencies['longitude'].value = 0.
      model.dependencies['latitude'].value = 0.
      model.dependencies['surface_downwelling_shortwave_flux'].value = 50.
      model.dependencies['surface_air_pressure'].value = 1.
      model.dependencies['wind_speed'].value = 5.
      model.dependencies['mole_fraction_of_carbon_dioxide_in_air'].value = 390.
      model.dependencies['number_of_days_since_start_of_the_year'].value = 1.
      
      # Verify the model is ready to be used
      model.cell_thickness=1.
      
      assert model.checkReady(), 'One or more model dependencies have not been fulfilled.'
      
      
      #variables with noise
#     varnames = ["Z5/c"]
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
      
      for i in range(len(t_eval)):
          dy = model.getRates()
          for j in range(len(model.state)):
              if i!=0:
                  if model.state_variables[j].name in varnames:
                      if i<15000:    #stop noise
                          y[s,i,j]=y[s,i-1,j]+dy[j]*dt  + D * sqrtdt * np.random.randn() * y[s,i-1,j]
                          A[s,i,j] = dy[j]
                      else:
                          y[s,i,j]=y[s,i-1,j]+dy[j]*dt
                          A[s,i,j] = dy[j]
                  else:
                          y[s,i,j]=y[s,i-1,j]+dy[j]*dt
                          A[s,i,j] = dy[j]
          model.state[:]=y[s,i,:]
if rank == 0:
    y_global = np.copy(y)
    for s in range(ntrials):
        r = s % nranks
        if r !=0 :
          rea = comm.recv(source=r, tag=11)
          y_global[int(rea['idx']),:,:] = rea['data']
          A[int(rea['idx']),:,:] = rea['grow']
        comm.Barrier()
else:
    for s in range(ntrials):
        r = s % nranks
        rea = {'idx':s, 'data': y[s,:,:], 'grow': A[s,:,:]}
        if r == rank :
            comm.send(rea, dest=0, tag=11)
        comm.Barrier()
comm.Barrier()

if rank == 0 :
  t_eval = t_eval / 86400
  y = np.copy(y_global)  
  


  
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
  n=np.arange(0,ntrials,1)
  sample_dim = f.createDimension('n', n.shape[0])
  
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
  
# temp = f.createVariable('temp', np.float64, ('time',))
# temp.units = 'C'
# temp.long_name = 'temperature in celsius'
# f.variables['temp'][:]=T

  sample = f.createVariable('n', np.float64, ('n',))
  sample.units = ''
  sample.long_name = 'sample'
  f.variables['n'][:]=n
  
  for v,variable in enumerate(model.state_variables):
     ncvar = variable.name.replace("/","_")
#    var = f.createVariable(ncvar, np.float64, ('time', 'z','lat','lon'))
     var = f.createVariable(ncvar, np.float64, ('n','time', 'z','lat','lon'))
     var.units = variable.units
     var.long_name = variable.long_name
#    f.variables[ncvar][:]=np.nanmean(y[:,:,v], axis=0)
     f.variables[ncvar][:]=y[:,:,v]


  # sensitivity targets
  laststeps = int(len(t_eval)/2.)
  y = np.nanmean(y,axis=0)
  kkk=0
  P1 = 0
  P1z = 1
  for v,variable in enumerate(model.state_variables):
      if variable.name.replace("/","_") == 'P1_c':
        P1 = v
  for v,variable in enumerate(model.state_variables):
     ncvar = variable.name.replace("/","_")
     if ncvar in varnames_:
         NA = (variation(y[-laststeps:,v])-variation(z[kkk,-laststeps:]))/D
         PN = np.mean( (y[-laststeps:,v]-np.mean(y[-laststeps:,v])) * (y[-laststeps:,P1]-np.mean(y[-laststeps:,P1])) ) / ( np.std(y[-laststeps:,v]) * np.std(y[-laststeps:,P1])  ) 
         cv = variation(y[-laststeps:,v])
         mean = np.mean(y[-laststeps:,v])
         if y[:,v].min() < 0.0001:
             ext=1
         else:
             ext=0
         
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
  
         ncfstate = variable.name.replace("/","_") + "_NA"
         P1_fstate = f.createVariable(ncfstate, np.float32, ('z','lat','lon'))
         P1_fstate.units = ''
         P1_fstate.long_name = variable.name.replace("/","_") + "_NA"
         f.variables[ncfstate][:]=NA
         print(ncfstate,NA)

         ncfstate = variable.name.replace("/","_") + "_PN"
         P1_fstate = f.createVariable(ncfstate, np.float32, ('z','lat','lon'))
         P1_fstate.units = ''
         P1_fstate.long_name = variable.name.replace("/","_") + "_PN"
         f.variables[ncfstate][:]=PN
         print(ncfstate,PN)
         print('*****')
  
         ncfstate = variable.name.replace("/","_") + "_ext"
         P1_fstate = f.createVariable(ncfstate, np.float32, ('z','lat','lon'))
         P1_fstate.units = ''
         P1_fstate.long_name = variable.name.replace("/","_") + "_ext"
         f.variables[ncfstate][:]=ext
  
  f.close()
