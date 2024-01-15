import numpy as np
import matplotlib.pyplot as plt





fig,ax = plt.subplots()

def logistic(a=4, n=100000, x0=0.4):
    x = np.zeros(n)
    x[0] = x0
    for i in range(n-1):
        x[i+1] = a*x[i]*(1-x[i])
    return(x)

interval = np.linspace(2.8,4.1,20)
print(interval)
time_series = [logistic(a=x) for x in interval]

for i in range(len(interval)):
        if interval[i] > 3:
            c = 'dodgerblue'
        if interval[i] > 3.56995:
            c = 'rebeccapurple'
        if interval[i] < 3:
            c = 'cyan'
        y = time_series[i][-1000:]
        for eq in y:
            ax.scatter(interval[i],eq,s=1,marker='o',facecolors='none', edgecolors=c)
            

ax.set_xlabel('a')
ax.set_ylabel('Logistic')

fig.savefig('logistic_map_bifurcation.png', format='png',dpi=350)




