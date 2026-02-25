import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import lyapunovV as lv

def lorenz(xyz, *, s=10, r=28, b=2.667):
    """
    Parameters
    ----------
    xyz : array-like, shape (3,)
       Point of interest in three-dimensional space.
    s, r, b : float
       Parameters defining the Lorenz attractor.

    Returns
    -------
    xyz_dot : array, shape (3,)
       Values of the Lorenz attractor's partial derivatives at *xyz*.
    """
    x, y, z = xyz
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return np.array([x_dot, y_dot, z_dot])


dt = 0.01
num_steps = 1000000

xyzs = np.empty((num_steps + 1, 3))  # Need one more for the initial values
xyzs[0] = (0., 1., 1.05)  # Set initial values
# Step through "time", calculating the partial derivatives at the current point
# and using them to estimate the next point
for i in range(num_steps):
    xyzs[i + 1] = xyzs[i] + lorenz(xyzs[i]) * dt

x = xyzs[:,0]

#compute lyapunov exponent on x timeseries

lyap = lv.LYAP(x)

##1- Wolf method:
#wolf = lyap.lyap_e()
#
#print(f"Wolf lyapunov: {wolf}")

#2- Palladin method
mean_tau,pall,mean_d0 = lyap.lyap_e_paladin(epsilon=1e-1,t_steps=1000)

print(f"Palladin lyapunov: {pall}")
print(f"Palladin mean tau= {mean_tau}, mean distance={mean_d0}")




