import pickle
import glob
import numpy as np
import netCDF4 as nc
from xml.dom import minidom
import matplotlib.pyplot as plt


files = ['small_cycles_result.nc','result_periodic.nc','result_chaos.nc']
#files = ['result_stationary.nc','result_periodic.nc','result_chaotic.nc']
conf = ['stationary','periodic','chaotic']
colors = ['cyan','dodgerblue','rebeccapurple']
amplitudes = [10,30,90,180,365]
files.reverse()
conf.reverse()
colors.reverse()
fig,axs = plt.subplots(3,3,sharex=True)#,constrained_layout=True)
plt.subplots_adjust(top=0.90)
plt.subplots_adjust(left=0.13)
plt.subplots_adjust(bottom=0.12)
axs = axs.ravel()
for inc,ncname in enumerate(files):
#   ncname = 'result.nc'
    f_noise = nc.Dataset(ncname)
    varnames=['B1_c','P1_c','P2_c','P3_c','P4_c','Z5_c','Z6_c','Z3_c','Z4_c']
    t = np.linspace(0, 7300., 100000)
    titles = ['Bacteria B1', 'Diatoms P1', 'Nanoflagellates P2', 'Picophytoplankton P3', 'Dinoflagellates P4', 'Microzooplankton Z5', 'het. Nanoflagellates Z6', 'carn. Mesozooplankton Z3', 'omn. Mesozooplankton Z4']
    for iax,ax in enumerate(axs):
    #   y = running_mean(f_noise.variables[varnames[iax]][:,0,0,0],410)
        ln = ax.plot(t[15000:],f_noise.variables[varnames[iax]][0,15000:,0,0,0],label=conf[inc],c=colors[inc])
    #   lns1 = ax.plot(t[-5000:],f_noise.variables[varnames[iax]][-5000:,0,0,0]-y[-5000:],label='noise',color='black')
        if inc == 0 and iax == 0:
            lns = ln 
        elif iax == 0:
            lns += ln
        ax.set_title(titles[iax],fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=8)

labels = [l.get_label() for l in lns]
fig.legend(lns, labels,ncol=3, loc='upper center',frameon=False)
fig.text(0.5, 0.04, 'time [days]', ha='center')
fig.text(0.04, 0.5, 'Biomass [$mg C/m^3$]', va='center', rotation='vertical')
fig.savefig('configurations.png', format='png',dpi=350)

fig1,axs1 = plt.subplots()
for inc,ncname in enumerate(files):
#   ncname = 'result.nc'
    f_noise = nc.Dataset(ncname)
    axs1.plot(f_noise.variables['sample'][:],f_noise.variables['max_lyapunov'][:,0,0,0],label=conf[inc],c=colors[inc])
axs1.axhline(y=0.001, color='k', linestyle='dotted')
axs1.set_yscale('log')
axs1.set_ylabel('Max Lyapunov Exp. [$d^{-1}$]')
#axs1.set_xlabel(r'$\tau1\;[days]$')
axs1.set_xlabel('Amplitude [$^\circ$C]')
axs1.set_xticks(f_noise.variables['sample'][::2])
fig1.legend(loc='center right', bbox_to_anchor=(0.4,0.8))
fig1.savefig('lyapunov.png', format='png',dpi=350)


#plot the timeseries with temperature fluctuations of 5C

fig,axs = plt.subplots(3,3,sharex=True)#,constrained_layout=True)
plt.subplots_adjust(top=0.90)
plt.subplots_adjust(left=0.13)
plt.subplots_adjust(bottom=0.12)
axs = axs.ravel()
for inc,ncname in enumerate(files):
#   ncname = 'result.nc'
    f_noise = nc.Dataset(ncname)
    varnames=['B1_c','P1_c','P2_c','P3_c','P4_c','Z5_c','Z6_c','Z3_c','Z4_c']
    t = np.linspace(0, 7300., 100000)
    titles = ['Bacteria B1', 'Diatoms P1', 'Nanoflagellates P2', 'Picophytoplankton P3', 'Dinoflagellates P4', 'Microzooplankton Z5', 'het. Nanoflagellates Z6', 'carn. Mesozooplankton Z3', 'omn. Mesozooplankton Z4']
    for iax,ax in enumerate(axs):
    #   y = running_mean(f_noise.variables[varnames[iax]][:,0,0,0],410)
        ln = ax.plot(t[15000:],f_noise.variables[varnames[iax]][5,15000:,0,0,0],label='$\lambda$= '+str(round(float(f_noise.variables['max_lyapunov'][5,0,0,0]),3)),c=colors[inc])
    #   lns1 = ax.plot(t[-5000:],f_noise.variables[varnames[iax]][-5000:,0,0,0]-y[-5000:],label='noise',color='black')
        if inc == 0 and iax == 0:
            lns = ln
        elif iax == 0:
            lns += ln
        ax.set_title(titles[iax],fontsize=8)
        ax.tick_params(axis='both', which='major', labelsize=8)

labels = [l.get_label() for l in lns]
fig.legend(lns, labels,ncol=3, loc='upper center',frameon=False)
fig.text(0.5, 0.04, 'time [days]', ha='center')
fig.text(0.04, 0.5, 'Biomass [$mg C/m^3$]', va='center', rotation='vertical')
fig.savefig('5C_configurations.png', format='png',dpi=350)

#interannaul variability of P3

fig,axs = plt.subplots(3,1,sharex=True)
axs = axs.ravel()
plt.subplots_adjust(top=0.90)
plt.subplots_adjust(left=0.14)
plt.subplots_adjust(right=0.80)
plt.subplots_adjust(bottom=0.12)
colors = ['b','g','r','c','m']
for inc,ncname in enumerate(files):
#   ncname = 'result.nc'
    f_noise = nc.Dataset(ncname)
    initial = 0
    for kk in range(int(len(f_noise.variables['Z6_c'][9,-30000:,0,0,0])/5000)-1):
        y = f_noise.variables['P1_c'][2,-25000:,0,0,0]
        time = np.linspace(0,365,5000)
        ln = axs[inc].plot(time,y[initial:initial+5000],label=str(kk),c=colors[kk])
        axs[inc].annotate(conf[inc], xy=(0.2,0.9), xycoords='axes fraction', xytext=(1, -1), textcoords='offset points', horizontalalignment='left', verticalalignment='top')
        if inc == 0 and kk == 0:
            lnp = ln
        elif inc == 0:
            lnp += ln
        if inc == 1:#'per' in conf[inc]:
        # inset axes....
#           x_detail = time[3972:4383]
            x_detail = time[3900:4200]
            y_detail = y[initial+3900:initial+4200]
            if kk == 0:
                sub_axes = plt.axes([.625,.45,.10,.10])
            sub_axes.plot(x_detail,y_detail,c=colors[kk])
#           plt.xticks(visible=False)
            plt.yticks(visible=False)
        if inc == 2:#'per' in conf[inc]:
        # inset axes....
#           x_detail = time[3972:4383]
            x_detail = time[3900:4200]
            y_detail = y[initial+3900:initial+4200]
            if kk == 0:
                sub_axes = plt.axes([.625,.20,.10,.10])
            sub_axes.plot(x_detail,y_detail,c=colors[kk])
#           plt.xticks(visible=False)
            plt.yticks(visible=False)
        initial += 5000
labels = [l.get_label() for l in lnp]
fig.legend(lnp, labels,title="year",loc='center left', bbox_to_anchor=(0.85,0.5))
fig.text(0.04, 0.5, 'Biomass P1 [$mg C/m^3$]', va='center', rotation='vertical')
fig.text(0.5, 0.04, 'time [days]', ha='center')
fig.savefig('annualvariability.png', format='png',dpi=350)


#plot temperature

T_tot = f_noise.variables['temp'][:]
MeanTemp = 15                       # Average temperature in the country
YearlyAmpl =  5    # Amplitude of the yearly cycle

# Total seconds in year
TotalHours = 24*365*60*60 #year period
tau        = 24*60*60     #day period

# Generate the frequency components of the data
t_eval = np.linspace(0, 7300.*86400, 100000)
YC = -YearlyAmpl*np.cos( (2*np.pi)*t_eval/TotalHours )

T = T_tot - MeanTemp - YC

fig,ax = plt.subplots()
#for kk in range(int(len(T[-30000:])/5000)-1):
#        y = T[-25000:]
#        time = np.linspace(0,365,5000)
#        ax.plot(time,T[initial:initial+5000],label=str(kk),c=colors[kk])
#        initial += 5000
ax.plot(t_eval/86400,T,c='k')
ax.set_xlabel('time [days]')
ax.set_ylabel('Amplitude [$^\circ$C]')
fig.savefig('noise_temperature.png', format='png',dpi=350)
