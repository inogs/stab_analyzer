import pickle
import glob
import numpy as np
import netCDF4 as nc
from xml.dom import minidom
import matplotlib.pyplot as plt


pkname = 'var.pickle'
infile = open(pkname,'rb')
new_dict = pickle.load(infile)
infile.close()
mean_var = new_dict.get('C')
var = new_dict.get('SC')
fig,axs = plt.subplots(3,3)#,constrained_layout=True)
fig.tight_layout()
plt.subplots_adjust(top=0.90)
plt.subplots_adjust(left=0.12)
plt.subplots_adjust(bottom=0.10)
plt.rcParams.update({'font.size': 5})
#fig.tight_layout()
axs = axs.ravel()
titles = ['Bacteria B1', 'Diatoms P1', 'Nanoflagellates P2', 'Picophytoplankton P3', 'Dinoflagellates P4', 'Microzooplankton Z5', 'het. Nanoflagellates Z6', 'carn. Mesozooplankton Z3', 'omn. Mesozooplankton Z4']
for iax,ax in enumerate(axs):
    hist,bin_edge = np.histogram(var[:,iax],bins=50)
    bins=np.zeros(len(hist))
    for b in range(len(bin_edge)-1):
        bins[b]=(bin_edge[b+1]+bin_edge[b])/2
    ax.scatter(bins, hist/len(var),s=0.1,c='black',label=f'not-steady mean variance = {round(mean_var[iax],3)}')
    ax.set_xlim([0.05,1])
    ax.set_ylim([0.0,0.01])
    ax.set_title(titles[iax],fontsize=5)
    ax.tick_params(axis='both', which='major', labelsize=5)
    string = 'mean = '+ str(round(mean_var[iax],3)).zfill(4)
    ax.annotate(string, xy=(0.6,0.9), xycoords='axes fraction', xytext=(1,-1), textcoords='offset points', horizontalalignment='right', verticalalignment='top')
fig.text(0.5, 0.04, 'Coefficient of Variation of Concentration ', ha='center', fontsize=8)
fig.text(0.04, 0.5, 'Normalized counts', va='center', rotation='vertical',fontsize=8)
#labels = [l.get_label() for l in ax]
print(np.sum(hist/len(var)))
#legend = [] 
#legend.append(ax.get_legend_handles_labels())
#legend.append(ax1.get_legend_handles_labels())

#fig.legend(lns, labels,ncol=2, loc='upper center')
#fig.tight_layout()
fig.savefig('FINAL_mean_var_cycles.png', format='png',dpi=150)

