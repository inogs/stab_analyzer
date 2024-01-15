import ordpy as od
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import pandas as pd

ncname = 'noise_result_stationary.nc'

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
      dt = t_eval[1]-t_eval[0]  # Time step.
      tt = t_eval[-1]  # Total time.
      n = int(tt / dt)+1  # Number of time steps.
      sigma_bis = 1/tau*np.sqrt(D)
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
#taux = np.arange(1,365,10)
taux = np.linspace(1,5000,100,dtype=int)
#for iv,var in enumerate(P1):
H = np.zeros(len(taux))
C = np.zeros(len(taux))
fig,axs = plt.subplots(1,2)
axs = axs.ravel()
for i in range(len(tau)):
    print(np.shape(T[i]))
    for it,tx in enumerate(taux):
        H[it],C[it] = od.complexity_entropy(T[i],taux=tx)
#   axs[2].plot(H,C,c=colors[i],label=str(tau[i]))
#   axs[3].plot(taux/5000*365,H,c=colors[i],linestyle='dashed',label=r'$H_{S}$')
#   axs[3].plot(taux/5000*365,C,c=colors[i],linestyle='dotted',label=r'$C_{JS}$')
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
        H[it],C[it] = od.complexity_entropy(P1i,taux=tx)
    axs[0].plot(H,C,c=colors[i],label=str(tau[i]))
    axs[1].plot(taux/5000*365,H,c=colors[i],linestyle='dashed',label=r'$H_{S}$')
    axs[1].plot(taux/5000*365,C,c=colors[i],linestyle='dotted',label=r'$C_{JS}$')
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
axs[1].set_ylabel(r'$H_{S},C_{JS}$')
axs[1].set_xlabel(r'$\epsilon$ $[days]$')
#axs[2].legend()
#axs[3].set_ylabel(r'$H_{S},C_{JS}$')
#axs[3].set_xlabel(r'embedding delay $[days]$')
#axs[3].legend()
axs[0].set_xlim(0.0,1.0)
#axs[1].set_xlim(0.0,1.0)
axs[0].set_ylim(0.0,0.5)
#axs[1].set_ylim(0.0,0.5)
fig.suptitle('P1 diatoms')
fig.tight_layout()
#plt.show()
fig.savefig('stochastic_plane.png',dpi=350)
