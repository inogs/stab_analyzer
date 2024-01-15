import ordpy as od
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd

files = ['result_chaos.nc','result_periodic.nc','noise_result_stationary.nc']
conf = ['endogenous chaos',r'induced chaos $5^\circ C$','stoch temp']
colors = ['k','b','r']

#HC for logistic map with parameter a=4 --> chaotic
H_logistic = 0.7387519146337446
C_logistic = 0.2918153116452781

fig,axs = plt.subplots(2,2,figsize=(5,6))
fig1,ax1 = plt.subplots(3,1)
axs = axs.ravel()
ax1 = ax1.ravel()
for ifc,ncname in enumerate(files):
        
    #ncname = 'noise_result_stationary.nc'
    
    f_noise = nc.Dataset(ncname)
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
    
    if 'noise' in ncname:
        P1 = T_tot[-80000:]
    else:
        P1 = f_noise.variables['P1_c'][5,-80000:,0,0,0]
    P1_1 = P1[:20000]
    P1_2 = P1[20000:40000]
    P1_3 = P1[40000:60000]
    P1_4 = P1[60000:]
    P1_1 = P1_2 - P1_1
    P1_2 = P1_3 - P1_4
    P1 = (P1_1+P1_2)/2
    ax1[ifc].plot(P1,c='k')
    taux = np.linspace(1,5000,100,dtype=int)
    H = np.zeros(len(taux))
    C = np.zeros(len(taux))
    for it,tx in enumerate(taux):
        H[it],C[it] = od.complexity_entropy(P1,taux=tx)
    axs[ifc].plot(H,C,c=colors[ifc],label=conf[ifc])
    axs[ifc].scatter(H_logistic,C_logistic,c='g')
    axs[3].plot(taux/5000*365,H,c=colors[ifc],linestyle='dashed',label=r'$H_{S}$')
    axs[3].plot(taux/5000*365,C,c=colors[ifc],linestyle='dotted',label=r'$C_{JS}$')
#   if 'periodic' in ncname:
#       P1 = f_noise.variables['P1_c'][19,-80000:,0,0,0]
#       P1_1 = P1[:20000]
#       P1_2 = P1[20000:40000]
#       P1_3 = P1[40000:60000]
#       P1_4 = P1[60000:]
#       P1_1 = P1_2 - P1_1
#       P1_2 = P1_3 - P1_4
#       P1 = (P1_1+P1_2)/2
#       taux = np.linspace(1,5000,100,dtype=int)
#       H = np.zeros(len(taux))
#       C = np.zeros(len(taux))
#       for it,tx in enumerate(taux):
#           H[it],C[it] = od.complexity_entropy(P1,taux=tx)
#       axs[0].plot(H,C,c='cyan',label=r'induced chaos $19^\circ C$')
#       axs[1].plot(taux/5000*365,H,c='cyan',linestyle='dashed',label=r'$H_{S}$')
#       axs[1].plot(taux/5000*365,C,c='cyan',linestyle='dotted',label=r'$C_{JS}$')


axs[0].set_ylabel(r'$C_{JS}$')
axs[1].set_ylabel(r'$C_{JS}$')
axs[2].set_ylabel(r'$C_{JS}$')
axs[0].set_xlabel(r'$H_{S}$')
axs[1].set_xlabel(r'$H_{S}$')
axs[2].set_xlabel(r'$H_{S}$')
axs[0].legend(loc='upper center', bbox_to_anchor=(0.5, 1.35),
          ncol=1)
axs[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.35),
          ncol=1)
axs[2].legend(loc='upper center', bbox_to_anchor=(0.5, 1.35),
          ncol=1)
axs[0].set_xlim(0.0,1.0)
axs[1].set_xlim(0.0,1.0)
axs[2].set_xlim(0.0,1.0)
axs[0].set_ylim(0.0,0.5)
axs[1].set_ylim(0.0,0.5)
axs[2].set_ylim(0.0,0.5)
axs[3].set_ylabel(r'$H_{S},C_{JS}$')
axs[3].set_xlabel(r'$\epsilon$ $[days]$')
axs[3].legend(loc='upper center', bbox_to_anchor=(0.5, 1.55),
          ncol=3)
fig.tight_layout()

fig.savefig('methods_chaos.png',dpi=350)
#plt.show()
