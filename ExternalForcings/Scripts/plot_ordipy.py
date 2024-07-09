import ordpy as od
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib

#ncname = 'noise_result_stationary.nc'
ncname = 'noise_result_chaotic.nc'

f_noise = nc.Dataset(ncname)
T_tot = f_noise.variables['temp'][:]
MeanTemp = 15                       # Average temperature in the country
YearlyAmpl =  5    # Amplitude of the yearly cycle

# Total seconds in year
TotalHours = 24*365*60*60 #year period
tau        = 24*60*60     #day period
tau_days=[10,30,90,180,365]
ntrials = len(tau_days)
# Generate the frequency components of the data
t_eval = np.linspace(0, 7300.*86400, 100000)
YC = -YearlyAmpl*np.cos( (2*np.pi)*t_eval/TotalHours )
T    = np.zeros((ntrials,len(YC)))
trials = 100
for s in range(ntrials):
      ##### add Ornstein-Uhlenbeck red noise to temperature
      sigma = 10.8    # Standard deviation.
      D     = 1.e5    # Noise amplitude [time]
      mu = 0.  # Mean
#     tau = tau_days*60*60/2/3.14 #Time constant
      tau = tau_days[s]*60*60/2/3.14 #Time constant
#     D = 2/tau/10**2
      dt = t_eval[1]-t_eval[0]  # Time step.
      tt = t_eval[-1]  # Total time.
      n = int(tt / dt)+1  # Number of time steps.
      sigma_bis = 1/tau*np.sqrt(D)
#     sigma_bis = np.sqrt(D/tau)
      sqrtdt = np.sqrt(dt)

      YearlyCycle =  np.zeros((ntrials,n))
      for i in range(ntrials):
        YearlyCycle[i,:] = YC
      F = np.zeros(trials)
      for i in range(n):
        # We update the process independently for
        # all trials
        F += dt * (-(F - mu) / tau) + \
            sigma_bis * sqrtdt * np.random.randn(trials)
        T[s,i]  = np.mean(F)  #temperature anomaly



T = T[:,-20000:]
P1 = f_noise.variables['P1_c'][:,-80000:,0,0,0]
tau = [10,30,90,180,365]
colors = list("royalblue cornflowerblue lightsteelblue lightsalmon tomato".split())
colors = ['#377eb8', '#ff7f00', '#4daf4a', '#f781bf', '#a65628']
#taux = np.arange(1,365,10)
taux = np.linspace(1,5000,100,dtype=int)
taux = np.logspace(np.log10(1),np.log10(5000),500,dtype=int)
#for iv,var in enumerate(P1):
H = np.zeros(len(taux))
C = np.zeros(len(taux))
fig,axs = plt.subplots(1,2)
axs = axs.ravel()
fig1,axs1 = plt.subplots(1,2)
axs1 = axs1.ravel()
fig3,ax3 = plt.subplots(1,2,figsize=(7,5))
fig4,ax4 = plt.subplots(1,3,figsize=(9,5),sharey=True)
for i in range(len(tau)):
    print(np.shape(T[i]))
    for it,tx in enumerate(taux):
        try:
            H[it],C[it] = od.complexity_entropy(T[i],taux=tx,dx=6)
        except:
            H[it],C[it] = np.nan, np.nan
#   axs1[0].plot(H,C,c=colors[i],label=str(tau[i]))
#   axs1[1].plot(taux/5000*365,H,c=colors[i],linestyle='dashed',label=r'$H_{S}$')
#   axs1[1].plot(taux/5000*365,C,c=colors[i],linestyle='dotted',label=r'$C_{JS}$')
    P1_1 = P1[i,:20000]
    P1_2 = P1[i,20000:40000]
    P1_3 = P1[i,40000:60000]
    P1_4 = P1[i,60000:]
    P1_1 = P1_2 - P1_1
    P1_2 = P1_3 - P1_4
    P1i = (P1_1+P1_2)/2
#   dic = {'signal':P1[i]}
#   data = pd.DataFrame(dic)
#   rolling_mean = data['signal'].rolling(window=410, center=True).mean()
#   P1i = data['signal'] - rolling_mean
    for it,tx in enumerate(taux):
        try:
            H[it],C[it] = od.complexity_entropy(P1i,taux=tx,dx=6)
        except:
            H[it],C[it] = np.nan, np.nan
    if i == 0:
        im0 = ax3[0].scatter(H,C,c=taux/5000*365,label=str(tau[i]),zorder=int(len(tau)-i),norm=matplotlib.colors.LogNorm(),cmap='Blues')
    if i == 1:
        im2 = ax4[0].scatter(H,C,c=taux/5000*365,label=str(tau[i]),zorder=int(len(tau)-i),norm=matplotlib.colors.LogNorm(),cmap='Oranges')
    if i == 2:
        im1 = ax3[1].scatter(H,C,c=taux/5000*365,label=str(tau[i]),zorder=int(len(tau)-i),norm=matplotlib.colors.LogNorm(),cmap='Greens')
    if i == 3:
        im3 = ax4[1].scatter(H,C,c=taux/5000*365,label=str(tau[i]),zorder=int(len(tau)-i),norm=matplotlib.colors.LogNorm(),cmap='Purples')
    if i == 4:
        im4 = ax4[2].scatter(H,C,c=taux/5000*365,label=str(tau[i]),zorder=int(len(tau)-i),norm=matplotlib.colors.LogNorm(),cmap='Reds')
    axs[0].plot(H[:-2],C[:-2],c=colors[i],label=str(tau[i]),zorder=int(len(tau)-i))
    axs[1].plot(taux[:-2]/5000*365,H[:-2],c=colors[i],linestyle='dashed',label=r'$H_{S}$')
    axs[1].plot(taux[:-2]/5000*365,C[:-2],c=colors[i],linestyle='dotted',label=r'$C_{JS}$')
    print(taux[:-10]/5000*365)

#plot max and min complexity
hmax, cmax = od.maximum_complexity_entropy(dx=6).T
hmin, cmin = od.minimum_complexity_entropy(dx=6, size=719).T

#fit to obtain max curve
import numpy.polynomial.polynomial as poly

coefs = poly.polyfit(hmax, cmax, 4)
ffit = poly.polyval(hmin, coefs)

hmax = hmin
cmax = ffit

axs[0].plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
axs[0].plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
axs[0].axvline(x=0.45,linestyle='--',c='k',linewidth=1.)
axs[0].axvline(x=0.70,linestyle='--',c='k',linewidth=1.)
axs[1].axhline(y=0.45,linestyle='--',c='k',linewidth=1.)
axs[1].axhline(y=0.70,linestyle='--',c='k',linewidth=1.)

axs[0].set_ylabel(r'$C_{JS}$')
axs[0].set_xlabel(r'$H_{S}$')
#axs[0].set_title(r'P1 diatoms')
axs[0].legend(title=r'$\tau$')
#axs[1].legend(title=r'$\tau$')
#axs[2].legend()
#axs[3].legend()
#axs[1].set_ylabel(r'$C_{JS}$')
#axs[1].set_xlabel(r'$H_{S}$')
#axs[1].set_title(r'Temperature Anomaly')
axs[1].set_ylabel(r'$C_{JS},H_{S}$')
axs[1].set_xlabel(r'$\epsilon$ $[days]$')
#axs[2].legend()
#axs[3].set_ylabel(r'$H_{S},C_{JS}$')
#axs[3].set_xlabel(r'embedding delay $[days]$')
#axs[3].legend()
axs[0].set_xlim(0.0,1.0)
#axs[1].set_xlim(0.0,1.0)
axs[0].set_ylim(0.0,0.5)
#axs[1].set_ylim(0.0,0.5)
axs[1].set_xscale('log')
fig.suptitle('P1 diatoms')
fig.tight_layout()
#plt.show()
fig.savefig('stochastic_plane.png',dpi=350)

fig,ax = plt.subplots()
ax.plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax.plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax.set_xlim(0.0,1.0)
ax.set_ylim(0.0,0.5)
ax.set_ylabel(r'$C_{JS}$')
ax.set_xlabel(r'$H_{S}$')
fig.savefig('empty_plane.pdf',dpi=500)

axs1[0].plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
axs1[0].plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
axs1[0].axvline(x=0.45,linestyle='--',c='k',linewidth=1.)
axs1[0].axvline(x=0.70,linestyle='--',c='k',linewidth=1.)
axs1[1].axhline(y=0.45,linestyle='--',c='k',linewidth=1.)
axs1[1].axhline(y=0.70,linestyle='--',c='k',linewidth=1.)

axs1[0].set_ylabel(r'$C_{JS}$')
axs1[0].set_xlabel(r'$H_{S}$')
axs1[0].legend(title=r'$\tau$')
axs1[1].set_ylabel(r'$C_{JS},H_{S}$')
axs1[1].set_xlabel(r'$\epsilon$ $[days]$')
axs1[0].set_xlim(0.0,1.0)
axs1[0].set_ylim(0.0,0.5)
axs1[1].set_xscale('log')
fig1.tight_layout()
#plt.show()
fig1.savefig('temp_stochastic_plane.png',dpi=350)

ax3[0].plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax3[0].plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax3[1].plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax3[1].plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax3[0].set_ylabel(r'$C_{JS}$')
ax3[0].set_xlabel(r'$H_{S}$')
ax3[1].set_xlabel(r'$H_{S}$')
ax3[0].set_title(r'$\tau=10\,days$')
ax3[1].set_title(r'$\tau=90\,days$')
#ax3.legend(title=r'$\tau$')
ax3[0].set_xlim(0.0,1.0)
ax3[0].set_ylim(0.0,0.5)
ax3[1].set_xlim(0.0,1.0)
ax3[1].set_ylim(0.0,0.5)
ax3[0].axvline(x=0.45,linestyle='--',c='k',linewidth=1.)
ax3[0].axvline(x=0.70,linestyle='--',c='k',linewidth=1.)
ax3[1].axvline(x=0.45,linestyle='--',c='k',linewidth=1.)
ax3[1].axvline(x=0.70,linestyle='--',c='k',linewidth=1.)
cbar = fig3.colorbar(im0, ax=ax3[0])
#cbar.set_label(r'$\epsilon$ $[days]$')
cbar1 = fig3.colorbar(im1, ax=ax3[1])
cbar1.set_label(r'$\epsilon$ $[days]$')
fig3.savefig('temp_stochastic_plane.png',dpi=500)

ax4[0].plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax4[0].plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax4[1].plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax4[1].plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax4[2].plot(hmin, cmin, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax4[2].plot(hmax, cmax, linewidth=1., color='#202020', zorder=0, alpha=0.4)
ax4[0].set_xlim(0.0,1.0)
ax4[0].set_ylim(0.0,0.5)
ax4[1].set_xlim(0.0,1.0)
ax4[1].set_ylim(0.0,0.5)
ax4[2].set_xlim(0.0,1.0)
ax4[2].set_ylim(0.0,0.5)
ax4[0].set_ylabel(r'$C_{JS}$')
ax4[0].set_xlabel(r'$H_{S}$')
ax4[1].set_xlabel(r'$H_{S}$')
ax4[2].set_xlabel(r'$H_{S}$')
ax4[0].set_title(r'$\tau=30\,days$')
ax4[1].set_title(r'$\tau=180\,days$')
ax4[2].set_title(r'$\tau=365\,days$')
ax4[0].axvline(x=0.45,linestyle='--',c='k',linewidth=1.)
ax4[0].axvline(x=0.70,linestyle='--',c='k',linewidth=1.)
ax4[1].axvline(x=0.45,linestyle='--',c='k',linewidth=1.)
ax4[1].axvline(x=0.70,linestyle='--',c='k',linewidth=1.)
ax4[2].axvline(x=0.45,linestyle='--',c='k',linewidth=1.)
ax4[2].axvline(x=0.70,linestyle='--',c='k',linewidth=1.)
cbar2 = fig4.colorbar(im2, ax=ax4[0])
#cbar2.set_label(r'$\epsilon$ $[days]$')
cbar3 = fig4.colorbar(im3, ax=ax4[1])
#cbar3.set_label(r'$\epsilon$ $[days]$')
cbar4 = fig4.colorbar(im4, ax=ax4[2])
cbar4.set_label(r'$\epsilon$ $[days]$')
fig4.savefig('add_temp_stochastic_plane.png',dpi=500)
