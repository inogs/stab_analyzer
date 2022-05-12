import pickle
import glob
import numpy as np
import netCDF4 as nc
from xml.dom import minidom
import matplotlib.pyplot as plt


pkname = 'std.pickle'
infile = open(pkname,'rb')
new_dict = pickle.load(infile)
infile.close()
cycles = new_dict.get('C')
linear = new_dict.get('NC')
t = np.linspace(0, 3650., 10000)
fig,axs = plt.subplots(3,3)#,constrained_layout=True)
fig.tight_layout()
plt.subplots_adjust(top=0.90)
plt.subplots_adjust(left=0.12)
plt.subplots_adjust(bottom=0.10)
#fig.tight_layout()
axs = axs.ravel()
titles = ['Bacteria B1', 'Diatoms P1', 'Nanoflagellates P2', 'Picophytoplankton P3', 'Dinoflagellates P4', 'Microzooplankton Z5', 'het. Nanoflagellates Z6', 'carn. Mesozooplankton Z3', 'omn. Mesozooplankton Z4']
for iax,ax in enumerate(axs):
    lns1 = ax.plot(t,cycles[iax],label='not-steady',color='black')
    ax1 = ax.twinx()
    lns2 = ax1.plot(t,linear[iax],label='steady',color='red')
    ax.set_title(titles[iax],fontsize=4)
    ax.tick_params(axis='both', which='major', labelsize=5)
    ax1.tick_params(axis='both', which='major', labelsize=5)
    ax1.tick_params(axis='y', colors='red')
    lns = lns1+lns2
    #labs = [l.get_label() for l in lns]
    #ax.legend(lns, labs, loc=0)
fig.text(0.5, 0.04, 'time [days]', ha='center')
fig.text(0.04, 0.5, 'Standard Deviation of the Concentration [$mg C/m^3$]', va='center', rotation='vertical')
labels = [l.get_label() for l in lns]
#legend = [] 
#legend.append(ax.get_legend_handles_labels())
#legend.append(ax1.get_legend_handles_labels())

fig.legend(lns, labels,ncol=2, loc='upper center')
#fig.tight_layout()
fig.savefig('FINAL_var_cycles.png', format='png',dpi=150)

shannon_cycles = new_dict.get('SC')/np.log(9)
shannon_lin = new_dict.get('SNC')/np.log(9)
print('shannon index cycles: ', np.mean(shannon_cycles))
print('shannon index no-cycles: ', np.mean(shannon_lin))

fig1,axs1 = plt.subplots()
hist,bin_edge = np.histogram(shannon_cycles)
bins=np.zeros(len(hist))
for b in range(len(bin_edge)-1):
    bins[b]=(bin_edge[b+1]+bin_edge[b])/2
axs1.scatter(bins, hist/len(shannon_cycles),c='black',label=f'not-steady mean Shannon index = {round(np.mean(shannon_cycles),3)}')

hist,bin_edge = np.histogram(shannon_lin)
bins=np.zeros(len(hist))
for b in range(len(bin_edge)-1):
    bins[b]=(bin_edge[b+1]+bin_edge[b])/2
axs1.scatter(bins, hist/len(shannon_lin),c='red',label=f'steady mean Shannon index = {round(np.mean(shannon_lin),3)}')
axs1.set_xlabel('normalized Shannon index')

fig1.legend()
fig1.savefig('FINAL_shannon.png', format='png',dpi=150)


