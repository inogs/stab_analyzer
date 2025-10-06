import lyapunovV as lv
import numpy as np
import matplotlib.pyplot as plt

def piecewise_map(ini,a,b,T):
    """ Piecewise linear map on [0,1].
    a > 1, 0 < b < 1, T integer     
    """
    if a <= 1 or b >= 1 or b <= 0 or T <= 0 or not isinstance(T, int):
        raise ValueError("Parameters must satisfy: a > 1, 0 < b < 1, T positive integer.")
    #time lenght bewteen 0 and 100*T with timestep of 1
    time = np.arange(0,100*T,1)
    #time = np.arange(0,100000,1)
    x = np.zeros(len(time))
    x[0] = ini
    for t in range(1,len(time)):
        # (2n-1)*T <= t < (2n)*T
        if ((t//T) % 2) == 0:
            x[t] = np.mod(a * x[t-1],1.0)
        else:
            x[t] = b*x[t-1]

    return x
T_array = [1,10,50,100,200,500,700,1000]
T_array = np.logspace(0.1,3,100,dtype=int)
T_array = np.unique(T_array)
ini = 0.1
a = 2
b = 2/3

fig,ax = plt.subplots(1,2)

for T in T_array:
#T = 10

    x = piecewise_map(ini,a,b,int(T))
    if T>20 and T < 60:
        ax[1].plot(x)
    
    lyap = lv.LYAP(x)

    #wolf
    try:
        lyap_wolf = lyap.lyap_e(dt=1)
        print("Wolf method: ", lyap_wolf)
    except:
        lyap_wolf = np.nan
        print("Wolf method failed")


    #vulpiani
    try:
        lyap_vulp = lyap.lyap_e_paladin(dt=1,delta0=1e-2, Delta=0.3)[-1]
        print("Vulpiani method: ", lyap_vulp)
    except:
        lyap_vulp = np.nan
        print("Vulpiani method failed")

    
    ax[0].scatter(T,lyap_wolf,c='k')
    ax[0].scatter(T,lyap_vulp,c='r')
    ax[0].set_xscale('log')

    fig.savefig('piecewise_map.png')