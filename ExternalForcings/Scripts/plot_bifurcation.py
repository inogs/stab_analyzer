import pickle
import glob
import numpy as np
import netCDF4 as nc
from xml.dom import minidom
import matplotlib.pyplot as plt
from numpy import array, sign, zeros
from scipy.interpolate import interp1d


def hl_envelopes_idx(s, dmin=1, dmax=1, split=False):
    """
    Input :
    s: 1d-array, data signal from which to extract high and low envelopes
    dmin, dmax: int, optional, size of chunks, use this if the size of the input signal is too big
    split: bool, optional, if True, split the signal in half along its mean, might help to generate the envelope in some cases
    Output :
    lmin,lmax : high/low envelope idx of input signal s
    """

    # locals min
    lmin = (np.diff(np.sign(np.diff(s))) > 0).nonzero()[0] + 1
#   lmin = (np.diff(np.sign(np.diff(-s))) < 0).nonzero()[0] + 1
    # locals max
    lmax = (np.diff(np.sign(np.diff(s))) < 0).nonzero()[0] + 1

    if split:
        # s_mid is zero if s centered around x-axis or more generally mean of signal
        s_mid = np.mean(s)
        # pre-sorting of locals min based on relative position with respect to s_mid
        lmin = lmin[s[lmin]<s_mid]
        # pre-sorting of local max based on relative position with respect to s_mid
        lmax = lmax[s[lmax]>s_mid]

    # global min of dmin-chunks of locals min
    lmin = lmin[[i+np.argmin(s[lmin[i:i+dmin]]) for i in range(0,len(lmin),dmin)]]
    # global max of dmax-chunks of locals max
    lmax = lmax[[i+np.argmax(s[lmax[i:i+dmax]]) for i in range(0,len(lmax),dmax)]]

    return lmin,lmax




#ncname = 'PF_result_chaos.nc'
ncname = 'result_chaos.nc'
#ncname = 'noise_result_periodic.nc'
#ncname = 'noise_result_stationary.nc'
amplitudes = np.linspace(0,20,20)
#amplitudes = [10,30,90,180,365]
f_noise = nc.Dataset(ncname)
varnames=['B1_c','P1_c','P2_c','P3_c','P4_c','Z5_c','Z6_c','Z3_c','Z4_c']
t = np.linspace(0, 7300., 100000)
fig,axs = plt.subplots(3,3,sharex=True)#,constrained_layout=True)
fig.tight_layout()
plt.subplots_adjust(top=0.90)
plt.subplots_adjust(left=0.13)
plt.subplots_adjust(bottom=0.12)
#fig.tight_layout()
axs = axs.ravel()
titles = ['Bacteria B1', 'Diatoms P1', 'Nanoflagellates P2', 'Picophytoplankton P3', 'Dinoflagellates P4', 'Microzooplankton Z5', 'het. Nanoflagellates Z6', 'carn. Mesozooplankton Z3', 'omn. Mesozooplankton Z4']
fig1,ax1 = plt.subplots(3,1,sharex=True)
fig2,ax2 = plt.subplots()
ax1 = ax1.ravel()
for iax,ax in enumerate(axs):  
#   if iax == 1:
#       break   
    print('doing ',varnames[iax])
    lyap = f_noise.variables['max_lyapunov'][:,0,0,0]
    for i in range(len(f_noise.variables[varnames[iax]][:15,0,0,0,0])):
        if lyap[i] > 0.001:
            c = 'rebeccapurple'
        else:
            c = 'dodgerblue'
        y = f_noise.variables[varnames[iax]][i,-50000:,0,0,0] #last three years timeseries
        lmin, lmax = hl_envelopes_idx(y)#,dmin=410)
        for upper in y[lmax]:
            ax.scatter(amplitudes[i],upper,s=1,marker='o',facecolors='none', edgecolors=c)
        for lower in y[lmin]:
            ax.scatter(amplitudes[i],lower,s=1,marker='o',facecolors='none', edgecolors=c)
        if varnames[iax] == 'P1_c':
            if i == 1 :
                ax1[0].plot(t[-50000:],y,c='k')
                ax1[0].set_title(r'$A_y=1.0^\circ C$')
                ax1[0].plot(t[-50000:][lmax],y[lmax],c='r')
                ax1[0].plot(t[-50000:][lmin],y[lmin],c='r')
            if i == 5 :
                ax1[1].plot(t[-50000:],y,c='k')
                ax1[1].set_title(r'$A_y=5.0^\circ C$')
                ax1[1].plot(t[-50000:][lmax],y[lmax],c='r')
                ax1[1].plot(t[-50000:][lmin],y[lmin],c='r')
            if i == 10 :
                ax1[2].plot(t[-50000:],y,c='k')
                ax1[2].set_title(r'$A_y=10.0^\circ C$')
                ax1[2].plot(t[-50000:][lmax],y[lmax],c='r')
                ax1[2].plot(t[-50000:][lmin],y[lmin],c='r')
            for upper in y[lmax]:
                    ax2.scatter(amplitudes[i],upper,s=1,marker='o',facecolors='none', edgecolors=c)
            for lower in y[lmin]:
                    ax2.scatter(amplitudes[i],lower,s=1,marker='o',facecolors='none', edgecolors=c)
    ax.set_title(titles[iax],fontsize=8)
    ax.set_xticks([5,10,15])
#   ax.set_xscale('log')
#   ax.tick_params(axis='both', which='major', labelsize=7)
axs[8].set_xticklabels([5,10,15])

#fig.text(0.5, 0.04, r'$\tau$ [$days$]', ha='center')
fig.text(0.5, 0.04, 'Amplitude [$^\circ$C]', ha='center')
fig.text(0.04, 0.5, 'Biomass [$mg C/m^3$]', va='center', rotation='vertical')
#legend = [] 
#legend.append(ax.get_legend_handles_labels())
#legend.append(ax1.get_legend_handles_labels())

#fig.tight_layout()
fig1.text(0.04, 0.5, 'P1 Biomass [$mg C/m^3$]', va='center', rotation='vertical')
ax1[2].set_xlabel('time [days]')
fig1.tight_layout(rect=(0.1, 0.1, 0.9, 0.9))


ax2.set_xlabel('Amplitude [$^\circ$C]')
ax2.set_ylabel('Biomass [$mg C/m^3$]')





fig.savefig('fluctuations.png', format='png',dpi=350)
fig1.savefig('methods_envelope.png', format='png',dpi=350)
fig2.savefig('methods_bifurcation.png', format='png',dpi=350)


plt.show()


