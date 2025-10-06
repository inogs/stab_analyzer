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
    time = np.arange(0,10*T,1)
    x = np.zeros(len(time))
    x[0] = ini
    for t in range(1,len(time)):
        # (2n-1)*T <= t < (2n)*T
        if ((t//T) % 2) == 0:
            x[t] = np.mod(a * x[t-1],1.0)
        else:
            x[t] = b*x[t-1]

    return x

ini = 0.1
a = 2
b = 2/3
T = 10

x = piecewise_map(ini,a,b,T)
plt.plot(x)
plt.savefig("piecewise_map.png")

lyap = lv.LYAP(x)

#wolf
try:
    lyap_wolf = lyap.lyap_e(dt=1)
    print("Wolf method: ", lyap_wolf)
except:
    print("Wolf method failed")


#vulpiani
try:
    lyap_vulp = lyap.lyap_e_paladin(dt=1,delta0=1e-2, Delta=0.3)
    print("Vulpiani method: ", lyap_vulp)
except:
    print("Vulpiani method failed")