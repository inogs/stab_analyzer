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

command='/g100_scratch/userexternal/gocchipi/BFM_5D/' + str(ensamble_counter).zfill(6) + '/'
#output of the model for all samples
ntrials = 4000
#ntrials = 4  
D = 0.001  # Standard deviation.

t_eval = np.linspace(0, 1825.*86400, 50000) 
y = np.zeros((ntrials,len(t_eval),54))
fig, ax = plt.subplots(3,3)
ax = ax.ravel()
fig1, ax1 = plt.subplots(3,3)
ax1 = ax1.ravel()
k=0

subprocess.run(["mkdir","-p",command])
subprocess.run(["cp","fabm.yaml",command])

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
      
      # Time derivative
      #def dy(t0, y):
      #    model.state[:] = y
      #    return model.getRates()
      
      #variables with noise
      varnames = ["B1/c","P1/c","P2/c","P3/c","P4/c","Z5/c","Z6/c","Z3/c","Z4/c"]
      
      # Time-integrate over 1000 days (note: FABM's internal time unit is seconds!)
      dt = (t_eval[-1]-t_eval[0])/len(t_eval)
      sqrtdt = np.sqrt(dt)
      y[s,0,:]=model.state[:]
      
      #print(model.state) 
      for i in range(len(t_eval)):
          dy = model.getRates()
          for j in range(len(model.state)):
              if i!=0:
                  if model.state_variables[j].name in varnames:
                      y[s,i,j]=y[s,i-1,j]+dy[j]*dt + D * sqrtdt * np.random.randn() * y[s,i-1,j]
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
  
  # We display the histogram for a few points in
   # time
  for i in range(len(t_eval)):
        for j in range(len(model.state)):
            if i in (0, 10000, 49000):# and s == ntrials-1:
              if model.state_variables[j].name in varnames:
                hist1, _ = np.histogram(y[:,i,j], bins=bins)
                ax[k].plot((_[1:] + _[:-1]) / 2, hist1/ntrials,
                #ax[k].plot((bins[1:] + bins[:-1]) / 2, hist1/ntrials,
                    {0: '-', 10000: '.', 49000: '-.', }[i],
                    label=f"t={i * dt / 86400:.2f}"+" d")
                ax[k].set_title(model.state_variables[j].name,fontsize=8)
                ax[8].legend()
                ax[k].tick_params(axis='both', which='major', labelsize=7)
                ax[k].tick_params(axis='both', which='minor', labelsize=7)
                ax[k].label_outer()
                #my_xticks = ax[k].get_xticks()
                #ax[k].set_xticks([my_xticks[0], my_xticks[-1]])#, visible=True, rotation="horizontal")
                k+=1
        k=0
  #get deterministic trajectory
  ncname = 'result_ode.nc'
  f = nc.Dataset(ncname)
  varnames_ = []
  for v in varnames:
   varnames_.append(v.replace("/","_"))
  z = np.zeros((len(varnames),len(t_eval)))
  for it,tname in enumerate(varnames_):  
            z[it] = f.variables[tname][:,0,0,0]
  for j in range(len(model.state)):
    if model.state_variables[j].name in varnames:
        ax1[k].plot(t_eval,np.mean(y[:,:,j],0),c='black',label='stoch')
        #ax1[k].plot(t_eval,np.amax(y[:,:,j],0),c='red')
        #ax1[k].plot(t_eval,np.amin(y[:,:,j],0),c='red')
        ax1[k].tick_params(axis='both', which='major', labelsize=7)
        ax1[k].tick_params(axis='both', which='minor', labelsize=7)
        ax1[k].plot(t_eval,z[k,:],c='blue',label='det')
        ax1[k].set_title(model.state_variables[j].name,fontsize=8)
        ax1[8].legend()
        k+=1
  fig.text(0.5, 0.01, 'Concentration [$mg C/m^3$]', ha='center')
  fig.text(0.01, 0.5, 'PDF', va='center', rotation='vertical')
  fig.savefig('stoch_histo.png', format='png',dpi=150)
  fig1.text(0.5, 0.01, 'time [days]', ha='center')
  #fig1.text(0.01, 0.5, 'Mean Concentration [$mg C/m^3$]', va='center', rotation='vertical')
  fig1.tight_layout()
  fig1.savefig('stoch_mean.png', format='png',dpi=150)
  #Relaxation time
  tau = np.zeros(len(varnames))
  for i in range(len(t_eval)):
      k=0
      for j in range(len(model.state)):
          if model.state_variables[j].name in varnames:
              tau[k]+=np.abs(np.mean(y[:,i,j],axis=0)-np.mean(y[:,-1,j],axis=0)) / np.abs(np.mean(y[:,0,j],axis=0)-np.mean(y[:,-1,j],axis=0))    
              k+=1
#  for j in range(len(model.state)):
#     if model.state_variables[j].name in varnames:
#        print(np.var(y[:,-1,j],axis=0))
  print(D)
  print(tau)
  print(y[:,0,50])
  ##t_eval = np.linspace(0, 3650.*86400, 300000) 
##sol = scipy.integrate.solve_ivp(dy, [0., 3650.*86400], model.state, t_eval=t_eval)
##y = scipy.integrate.odeint(dy, model.state, t*86400)
#t = t_eval/86400
#
## Plot results
##import pylab
##pylab.plot(t, y)
##pylab.legend([variable.path for variable in model.state_variables])
##pylab.show()
#
#
#Nt=t.shape[0]
#deltaT=t[1]-t[0]
#laststeps = int(Nt/10) #compute the indicators just for this last timesteps
#freq = (1/deltaT) * np.linspace(0,laststeps/2,int(laststeps/2)) / laststeps
#
# Save results


fileoutput = 'result.nc'
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
## Sensibility targets
#
##critical lenght for cycles
#epsilon = 0.001
#for v,variable in enumerate(model.state_variables):
#
##coeff of variance
#    if math.isinf(variation(y[-laststeps:,v])) : 
#        p1_var=int(1)
#    elif variation(y[-laststeps:,v]) > epsilon :
#        p1_var = int(1)
#    else :
#        p1_var = int(0)
#
##Relative peaks in the fourier spac
#    PF = fft(y[-laststeps:,v])
#
#    PF2 = np.abs(PF / laststeps) #full spectrum  /Npoints to get the true amplitudes
#    PF1 = PF2[:laststeps//2]              #half spectrum
#    PF1[1:] = 2*PF1[1:]           #actual amplitude
#
#    #find peaks
#    peak_index0, dict_vals = find_peaks(np.concatenate(([min(PF1)],PF1,[min(PF1)])), distance=5)
#    #find the biggest two
#    peak_index0 = peak_index0-1
#    heights0 = np.sort(PF1[peak_index0])[::-1]
#    if len(heights0) >1 :
#        if heights0[0] >0.001 :
#            p1_peak = heights0[1]/heights0[0]
#        else :
#            p1_peak = 0
#    else :
#        p1_peak = 0
#    #make it discrete
#    if p1_peak > epsilon :
#        p1_peak = int(1)
#    else :
#        p1_peak = int(0)
#
##fourier spectrum plot
##    if v==1:
##        pylab.plot(freq,PF1)
##        pylab.show()
#
##Relaxation time
##    final_times = int( Nt/10)
##    yfinP = np.mean(y[-final_times:,v])
##    for i,Y in enumerate(y[-laststeps:,v]):
##        p1_tau =+ np.abs( Y - yfinP ) / np.abs(y[0,v]-yfinP)
#
##Autocorrelation
#
##    p1_acf = sm.tsa.acf(y[-laststeps:,v], nlags=1)[1]
##Last state
##    fstate = y[-1,v]
#
#    ncvar = variable.name.replace("/","_") + "_var"
#    P1_var = f.createVariable(ncvar, np.float32, ('z','lat','lon'))
#    P1_var.units = ' '
#    P1_var.long_name = ncvar + 'iance'
#    f.variables[ncvar][:]=p1_var
#
#    ncpeak = variable.name.replace("/","_") + "_peak"
#    P1_peak = f.createVariable(ncpeak, np.float32, ('z','lat','lon'))
#    P1_peak.units = ' '
#    P1_peak.long_name = ncpeak +'_height'
#    f.variables[ncpeak][:]=p1_peak
#
##    nctau = variable.name.replace("/","_") + "_tau"
##    P1_tau = f.createVariable(nctau, np.float32, ('z','lat','lon'))
##    P1_tau.units = 'days'
##    P1_tau.long_name = variable.name.replace("/","_") + "_relaxation_time"
##    f.variables[nctau][:]=p1_tau
#
##    nctau = variable.name.replace("/","_") + "_acf"
##    P1_acf = f.createVariable(nctau, np.float32, ('z','lat','lon'))
##    P1_acf.units = ''
##    P1_acf.long_name = variable.name.replace("/","_") + "_lag1_autocorrelation"
##    f.variables[nctau][:]=p1_acf
#    
##    ncfstate = variable.name.replace("/","_") + "_fstate"
##    P1_fstate = f.createVariable(ncfstate, np.float32, ('z','lat','lon'))
##    P1_fstate.units = ''
##    P1_fstate.long_name = variable.name.replace("/","_") + "_final_state"
##    f.variables[ncfstate][:]=fstate
f.close()
#subprocess.run(["cp","result.nc",command])
