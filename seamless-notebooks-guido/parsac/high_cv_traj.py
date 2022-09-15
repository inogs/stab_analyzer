import pickle
import glob
import numpy as np
import netCDF4 as nc
from xml.dom import minidom
from scipy.stats import variation 
try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False

print("Init ...", flush=True)
print(isParallel, flush=True)
#get all nc paths
filenames = glob.glob('/g100_scratch/userexternal/gocchipi/BFM_TOTAL/*/*.nc')
filenames.sort(key=lambda x: int(''.join(filter(str.isdigit, x)))) #sort by number

indexes = np.zeros(len(filenames))
for iv,var in enumerate(filenames):
     indexes[iv] = int(var.rsplit('/')[-2])

#interested variables names
varnames = ["B1_c","P1_c","P2_c","P3_c","P4_c","Z5_c","Z6_c","Z3_c","Z4_c"]
#get the cycling simulation for small N1p

pknamein = "bfm_sensitivity_total.pickle"
pknameout = 'bfm_sensitivity_total_script.pickle'
infile = open(pknamein,'rb')
new_dictin = pickle.load(infile)
infile.close()
inputs = new_dictin.get('X')
infile = open(pknameout,'rb')
new_dictout = pickle.load(infile)
infile.close()
outputs = new_dictout.get('Y')
outputs = outputs.tolist()
#restrict inputs to the existing directories
inputs = [inputs[int(i)] for i in indexes]
#get target names

mydoc = minidom.parse('bfm_sensitivity.xml')
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)

#remove 05 03c 03h R3c
keyo2 = 'O2_o'
key = 'O5'
keyOc= 'O3_c'
keyOh= 'O3h_h'
keyR3c = 'R3_c'
false_index = []
for i,name in enumerate(target):
    if name.find(keyo2) !=-1 or name.find(key) !=-1 or name.find(keyOc) !=-1 or name.find(keyOh) !=-1 or name.find(keyR3c) !=-1:
        false_index.append(i)
for i in range(len(outputs)):  #remove O5
    for ind in sorted(false_index, reverse=True):
         outputs[i].pop(ind)
for ind in sorted(false_index, reverse=True):
    target.pop(ind)
#get parameters name
#get parameters name from the xml
mydoc = minidom.parse('bfm_sensitivity.xml')
items = mydoc.getElementsByTagName('parameter')
parameters = []
for i in range(len(items)):
    string = items[i].attributes['variable'].value
    if i<200:
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
    else :
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
#get N1p index
N1p = 0
for i,name in enumerate(parameters):
    if name == "N1p":
        N1p = i

#count cycles        
count = np.zeros(len(inputs))
for o in range(len(outputs)):
    for i in range(int(len(outputs[o])/2)):
        if outputs[o][2*i]+outputs[o][2*i+1]==2:
            count[o]=1
        else :
            pass
#select directories with cycling behaviour for small N1p (<0.35)
cyclesind =[]
for i in range(len(inputs)):
    if count[i]==1 and inputs[i][N1p]<0.35 : 
        cyclesind.append(i)
#select directories with non-cycling behaviour for large N1p (>0.35)
linind = []
for i in range(len(inputs)):
    if count[i]==0 and inputs[i][N1p]>0.35 :
        linind.append(i)
#lenght of trajectories
lenght = 10000

#restirct to large covariance trajectories
pkname = 'var.pickle'
infile = open(pkname,'rb')
new_dict = pickle.load(infile)
infile.close()
cv = new_dict.get('SC')
idx = []
for iv in range(len(cv)):
        if np.max(cv[iv]) > 0.05: #collect trajectories with CV bigger than 0.05
                idx.append(iv)
print(idx,flush=True)

# cycling trajectories
output_cycles = np.zeros((len(varnames),lenght))
filenames_cycle = [filenames[i] for i in cyclesind]
filenames_cycle_var = [filenames_cycle[i] for i in idx]

p_cycles = np.zeros(len(varnames))
var_cycles = np.zeros((len(filenames_cycle_var),len(varnames))) 
count_s_c=0
sh=0
tot = np.zeros(len(varnames))
traj = np.zeros((len(filenames_cycle_var),len(varnames),lenght)) 
print(len(filenames_cycle_var))
for inc,ncname in enumerate(filenames_cycle_var) :
    i= inc % nranks
    if i == rank :
        f = nc.Dataset(ncname)
        for it,tname in enumerate(varnames):    #keep same order of xml file
            var = f.variables[tname][:,:,:,:]
            traj[inc,it,:] = var[:,0,0,0]
print('collecting ranks', flush=True)
if rank == 0:
    for ipc,ncname in enumerate(filenames_cycle_var):
        i = ipc % nranks
        if i!=0: 
            dic = comm.recv( source=i, tag=4 )
            traj_global[int(dic['idx']),:] = dic['data']
else :
        print('sending...',flush=True)
        for ipc,ncname in enumerate(filenames_cycle_var) :
            i= ipc % nranks
            dic = {'idx':ipc, 'data': traj[ipc,:,:]}
            if i== rank :
                comm.send( dic, dest=0, tag=4)
print('End of the loop',flush=True)
#compute means
if rank == 0 :
    print('writing pickle', flush=True)
    norm = 1.0 / len(filenames_cycle_var)
    pkname = 'traj.pickle'
    outfile = open(pkname,'wb')
    out_dict = {'T': traj_global}
    pickle.dump(out_dict,outfile)
    outfile.close()





