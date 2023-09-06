import pickle
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import argparse
from xml.dom import minidom
import glob

xml = '/g100_work/OGS21_PRACE_P/gocchipinti/PeriodicForcing/seamless-notebooks-guido/parsac/bfm_N1p.xml'
pkl = '/g100_work/OGS21_PRACE_P/gocchipinti/PeriodicForcing/seamless-notebooks-guido/parsac/bfm_NF.pickle'

#get parameters name from the xml
mydoc = minidom.parse(xml)
items = mydoc.getElementsByTagName('parameter')
parameters = [] 
for i in range(len(items)):
    string = items[i].attributes['variable'].value
    parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])

#get targets names
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)

infile = open(pkl,'rb')
new_dict = pickle.load(infile)
infile.close()
inputs = new_dict.get('X')
outputs = new_dict.get('Y')
#outputs = outputs.tolist()

threshold = [1.e-6,1.e-5,1.e-4,1.e-3,1.e-2,1.e-1]
lll = np.zeros(len(threshold))
for ie,epsilon in enumerate(threshold):
    count = np.zeros(len(inputs))
    #epsilon = 0.001
    for o in range(len(outputs)):
        for i in range(len(outputs[o])-2):
            if outputs[o][i] > epsilon:
                count[o]=1
            else :
                pass
#   lll = 0
    for l in count:
        if l==1:
            lll[ie] +=1
count = np.zeros(len(inputs))
epsilon = 0.001
for o in range(len(outputs)):
    for i in range(len(outputs[o])-2):
        if outputs[o][i] > epsilon:
            count[o]=1
        else :
            pass
lll = lll/len(outputs) *100
print('total sample',len(outputs))
print('non-stationary',lll)
print('fraction',lll/len(outputs))
fig, axs = plt.subplots(1, len(inputs[0]),sharey=True)
axs = axs.ravel()
#plt.rcParams.update({'font.size': 4})
x = np.zeros(len(inputs))
for i in range(len(inputs[0])):
    for j in range(len(inputs)):
        x[j]=inputs[j][i]
    hist,bin_edge = np.histogram(x, weights=count)
    hist_norm,bin_ed = np.histogram(x)
    bins=np.zeros(len(hist))
    for b in range(len(bin_edge)-1):
        bins[b]=(bin_edge[b+1]+bin_edge[b])/2
    axs[i].plot(bins, hist/hist_norm, c='black')

for iax,ax in enumerate(axs):
    ax.label_outer()
    ax.axhline(y=0.5, color='r', linestyle="dotted")
    ax.set_ylim(0, 1)
#   ax.set_xticks([])
#   ax.set_yticks([])
    ax.annotate(parameters[iax], xy=(0,1), xycoords='axes fraction', xytext=(1, -1), textcoords='offset points', horizontalalignment='left', verticalalignment='top')


fig.savefig('nonstat_frequency.png', format='png',dpi=150)

fig,ax = plt.subplots()
ax.plot(threshold,lll,c='k',ls='-',marker='2', markerfacecolor='none')
ax.set_xscale("log")
ax.set_xlabel(r'Threshold $\epsilon$')
ax.set_ylabel(r'Non-stationary frequency $[\%]$')
fig.savefig('threshold.png', format='png',dpi=150)


#lyapunov

lyap = outputs[:,-1]
ext  = outputs[:,-2]
cv_max = np.max(outputs[:,:-2],axis=1)
id_chaos = []
id_periodic = []
id_stationary = []
id_extincted = []
for i in range(len(lyap)):
    if ext[i] == 1:
        id_extincted.append(i)
    elif cv_max[i] < 1.e-3:
        id_stationary.append(i)
    elif lyap[i] < 1.e-3:
        id_periodic.append(i)
    else:
        id_chaos.append(i)

fig, ax = plt.subplots()

chosen_points_x = [0.06,0.08,0.15]
chosen_points_y = [0.6,0.6,0.6]

ax.axhline(y=0.6, color='k', linestyle=(0, (5, 10)),alpha=0.3)
ax.scatter(inputs[id_extincted,0],inputs[id_extincted,1],c='black', label='extincted',s=4)
ax.scatter(inputs[id_stationary,0],inputs[id_stationary,1],c='cyan', label='stationary',s=4)
ax.scatter(inputs[id_periodic,0],inputs[id_periodic,1],c='dodgerblue', label='periodic',s=4)
ax.scatter(inputs[id_chaos,0],inputs[id_chaos,1],c='darkorchid', label='chaotic',s=4)

ax.scatter(chosen_points_x,chosen_points_y,marker='*',c='magenta',s=28)

ax.set_xlabel('$PO_4$')
ax.set_ylabel(r'$\beta_z$')
ax.set_aspect('equal')
#fig.colorbar(h)
fig.legend(loc='center left', bbox_to_anchor=(0.7,0.5))
plt.show()
fig.savefig('scatter.png', format='png',dpi=350)

#rate of extinction

print('rate of extinction: ', len(id_extincted)/len(lyap)*100)

