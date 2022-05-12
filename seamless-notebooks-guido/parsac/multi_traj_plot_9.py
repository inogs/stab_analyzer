import pickle
import glob
import numpy as np
import netCDF4 as nc
from xml.dom import minidom
import matplotlib.pyplot as plt


pkname = 'traj.pickle'
infile = open(pkname,'rb')
new_dict = pickle.load(infile)
infile.close()
traj = new_dict.get('T')

#compute shannon index
sh = []
tot = np.zeros(len(traj[0]))
lenght = len(traj[0,0])
p_cycles = np.zeros((len(traj),len(traj[0])))
for i in range(len(traj)):
    for j in range(len(traj[0])):
        if np.mean(traj[i,j,-int(lenght/10):])<0:
            p_cycles[i,j] = 0
        else :
            p_cycles[i,j] = np.mean(traj[i,j,-int(lenght/10):])
for i in range(len(p_cycles)):
    for ip,pv in enumerate(p_cycles[i]):
        if pv ==0:
            tot[ip]=0
        else :
            tot[ip] = pv/np.sum(p_cycles[i])*np.log(pv/np.sum(p_cycles[i]))
    sh.append(-np.sum(tot))
print(sh/np.log(9))
print(np.shape(p_cycles))
print(np.mean(sh))
print(np.mean(sh)/np.log(9))

t = np.linspace(0, 3650., 10000)
for idx in range(len(traj)):
    fig,axs = plt.subplots(3,3)#,constrained_layout=True)
    fig.tight_layout()
    plt.subplots_adjust(top=0.90)
    plt.subplots_adjust(left=0.12)
    plt.subplots_adjust(bottom=0.10)
    #fig.tight_layout()
    axs = axs.ravel()
    titles = ['Bacteria B1', 'Diatoms P1', 'Nanoflagellates P2', 'Picophytoplankton P3', 'Dinoflagellates P4', 'Microzooplankton Z5', 'het. Nanoflagellates Z6', 'carn. Mesozooplankton Z3', 'omn. Mesozooplankton Z4']
    for iax,ax in enumerate(axs):
        lns1 = ax.plot(t,traj[idx,iax,:],label='not-steady',color='black')
        ax.set_title(titles[iax],fontsize=4)
        ax.tick_params(axis='both', which='major', labelsize=5)
    fig.text(0.5, 0.04, 'time [days]', ha='center')
    fig.text(0.04, 0.5, 'Mean Value Concentration [$mg C/m^3$]', va='center', rotation='vertical')
    stri = str(round(sh[idx]/np.log(9),3)).zfill(3)
    stri = 'shannon: '+stri
    plt.suptitle(stri,fontsize=7)
    string = 'traj_cycles_'+str(idx)+'.png'
    
    fig.savefig(string, format='png',dpi=150)

fig1,ax1 = plt.subplots()
for iv,var in enumerate(traj[2]):
    ax1.plot(t,var,label=titles[iv])
ax1.axis('off')
fig1.savefig('header.png', format='png',dpi=150)


pkname = 'means.pickle'
infile = open(pkname,'rb')
new_dict = pickle.load(infile)
infile.close()
linear = new_dict.get('NC')

mean_cycles=np.zeros(np.shape(traj[0]))
for iv,var in enumerate(traj):
    for i in range(len(var)):
        mean_cycles[i] += var[i]

mean_cycles = mean_cycles /len(traj)


fig2,axs = plt.subplots(3,3)#,constrained_layout=True)
fig2.tight_layout()
plt.subplots_adjust(top=0.90)
plt.subplots_adjust(left=0.12)
plt.subplots_adjust(bottom=0.10)
axs = axs.ravel()
for iax,ax in enumerate(axs):
  lns1 = ax.plot(t,mean_cycles[iax],label='not-steady',color='black')
  ax1 = ax.twinx()
  lns2 = ax1.plot(t,linear[iax],label='steady',color='red')
  ax.set_title(titles[iax],fontsize=4)
  ax.tick_params(axis='both', which='major', labelsize=5)
  ax1.tick_params(axis='y', colors='red')
  ax1.tick_params(axis='both', which='major', labelsize=5)
  lns = lns1+lns2
fig.text(0.5, 0.04, 'time [days]', ha='center')
fig.text(0.04, 0.5, 'Mean Value Concentration [$mg C/m^3$]', va='center', rotation='vertical')
labels = [l.get_label() for l in lns]
fig2.legend(lns, labels,ncol=2, loc='upper center')
#fig.tight_layout()
fig2.savefig('FINAL_mean_largefluct.png', format='png',dpi=150)
