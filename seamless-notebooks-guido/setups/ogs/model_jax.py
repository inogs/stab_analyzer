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
import jax

ensamble_counter=int(np.loadtxt('ensamble_counter.txt'))

command='/g100_scratch/userexternal/gocchipi/BFM_TOTAL/' + str(ensamble_counter).zfill(6) + '/'
subprocess.run(["mkdir","-p",command])
subprocess.run(["cp","fabm.yaml",command])

# Create model (loads fabm.yaml)
model = pyfabm.Model('fabm.yaml_cycles')

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
model.dependencies['wind_speed'].value = 1.
model.dependencies['mole_fraction_of_carbon_dioxide_in_air'].value = 390.
model.dependencies['number_of_days_since_start_of_the_year'].value = 1.

# Verify the model is ready to be used
model.cell_thickness=1.

assert model.checkReady(), 'One or more model dependencies have not been fulfilled.'

# Time derivative
def dy(y, t0):
    model.state[:] = y
    return model.getRates()

N = len(model.state)
# Compute Jacobian
def dyj(y):
    rates = jax.numpy.zeros(N)
    for i in range(N):
        rates = rates.at[i].set(model.getRates()[i])
    return rates 
def fL(SL,t):
#    res=jax.numpy.dot(jax.jacfwd(dyj)(y[i,:]),SL)
    f,res=jax.jvp(dyj,(y[i,:],),(y[i,:],))
    return res
states = jax.numpy.zeros(N)
for i in range(N):
    states = states.at[i].set(model.state[i])
T = np.linspace(0, 3650.*86400, 300000)#jax.numpy.zeros(NSTEPS)
NSTEP = len(T)
y = jax.numpy.zeros((NSTEP,N))
yL = jax.numpy.zeros((N,N))
stateL = jax.numpy.zeros((N,N))
ZNORM = jax.numpy.zeros(N)
GSC = jax.numpy.zeros(N)
CUM = jax.numpy.zeros((N,NSTEP))
CUM_NORM = jax.numpy.zeros((N,NSTEP))
#initial conditions for non linear system
y = y.at[0,:].set(model.state)
#initial conditions for linear system
for i in range(N):
    yL = yL.at[i,i].set(1.0)
#integration
for i in range(len(T)):
    XEND = T[i+1]
    t = [T[i],T[i+1]]
    for j in range(N):
        stateL = stateL.at[j,:].set(scipy.integrate.odeint(fL, yL[0,:], t)[-1])
    yL = jax.numpy.copy(stateL)
    y = y.at[i,:].set(scipy.integrate.odeint(dy, model.state, t)[-1])
    #normalize vectors
    for k in range(N):
    #    if k==0:
    #        ZNORM=ZNORM.at[k].set(jax.numpy.linalg.norm(yL[k,:]))
    #        yL=yL.at[k,:].set(yL[k,:]/ZNORM[k])
    #    else:
            a = yL[k]
            for j in range(N):
                if j<k:
                    a = a - jax.numpy.dot(yL[k,:],yL[j,:])*yL[j,:]
            print(a)
            yL = yL.at[k].set(a)
            ZNORM=ZNORM.at[k].set(jax.numpy.linalg.norm(yL[k,:]))
            yL=yL.at[k,:].set(yL[k,:]/ZNORM[k])
    #update running vector magnitudes
    if i==0:
        for k in range(N):
            CUM = CUM.at[k,i].set(np.log(ZNORM[k]))
    else:
        for k in range(N):
            CUM = CUM.at[k,i].set(CUM[k,i-1] + np.log(ZNORM[k]))
    #normalize exponent
    for k in range(N):
        CUM_NORM = CUM_NORM.at[k,i].set(CUM[k,i]/XEND)

#    print(CUM_NORM[:,i])            
        
t = T/86400
# Time-integrate over 1000 days (note: FABM's internal time unit is seconds!)
#t = np.linspace(0, 3650., 300000)
#y = scipy.integrate.odeint(dy, model.state, t*86400)

# Plot results
#import pylab
#pylab.plot(t, y)
#pylab.legend([variable.path for variable in model.state_variables])
#pylab.show()


Nt=t.shape[0]
deltaT=t[1]-t[0]
laststeps = int(Nt/10) #compute the indicators just for this last timesteps
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

for v,variable in enumerate(model.state_variables):
   ncvar = variable.name.replace("/","_")
   var = f.createVariable(ncvar, np.float64, ('time', 'z','lat','lon'))
   var.units = variable.units
   var.long_name = variable.long_name
   f.variables[ncvar][:]=y[:,v]

# Sensibility targets

#critical lenght for cycles
epsilon = 0.001
for v,variable in enumerate(model.state_variables):

#coeff of variance
    if math.isinf(variation(y[-laststeps:,v])) : 
        p1_var=int(1)
    elif variation(y[-laststeps:,v]) > epsilon :
        p1_var = int(1)
    else :
        p1_var = int(0)

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

#fourier spectrum plot
#    if v==1:
#        pylab.plot(freq,PF1)
#        pylab.show()

#Relaxation time
#    final_times = int( Nt/10)
#    yfinP = np.mean(y[-final_times:,v])
#    for i,Y in enumerate(y[-laststeps:,v]):
#        p1_tau =+ np.abs( Y - yfinP ) / np.abs(y[0,v]-yfinP)

#Autocorrelation

#    p1_acf = sm.tsa.acf(y[-laststeps:,v], nlags=1)[1]
#Last state
#    fstate = y[-1,v]

    ncvar = variable.name.replace("/","_") + "_var"
    P1_var = f.createVariable(ncvar, np.float32, ('z','lat','lon'))
    P1_var.units = ' '
    P1_var.long_name = ncvar + 'iance'
    f.variables[ncvar][:]=p1_var

    ncpeak = variable.name.replace("/","_") + "_peak"
    P1_peak = f.createVariable(ncpeak, np.float32, ('z','lat','lon'))
    P1_peak.units = ' '
    P1_peak.long_name = ncpeak +'_height'
    f.variables[ncpeak][:]=p1_peak

    nclyap = variable.name.replace("/","_") + "_lyap"
    P1_lyap = f.createVariable(nclyap, np.float32, ('z','lat','lon'))
    P1_lyap.units = 's-1'
    P1_lyap.long_name = ncpeak +'_lyapunov'
    f.variables[ncpeak][:]=CUM_NORM[v,-1]

#    nctau = variable.name.replace("/","_") + "_tau"
#    P1_tau = f.createVariable(nctau, np.float32, ('z','lat','lon'))
#    P1_tau.units = 'days'
#    P1_tau.long_name = variable.name.replace("/","_") + "_relaxation_time"
#    f.variables[nctau][:]=p1_tau

#    nctau = variable.name.replace("/","_") + "_acf"
#    P1_acf = f.createVariable(nctau, np.float32, ('z','lat','lon'))
#    P1_acf.units = ''
#    P1_acf.long_name = variable.name.replace("/","_") + "_lag1_autocorrelation"
#    f.variables[nctau][:]=p1_acf
    
#    ncfstate = variable.name.replace("/","_") + "_fstate"
#    P1_fstate = f.createVariable(ncfstate, np.float32, ('z','lat','lon'))
#    P1_fstate.units = ''
#    P1_fstate.long_name = variable.name.replace("/","_") + "_final_state"
#    f.variables[ncfstate][:]=fstate
f.close()
subprocess.run(["cp","result.nc",command])
