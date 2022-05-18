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
# cycling trajectories
output_cycles = np.zeros((len(varnames),lenght))
filenames_cycle = [filenames[i] for i in cyclesind]
p_cycles = np.zeros(len(varnames))
var_cycles = np.zeros((len(filenames_cycle),len(varnames))) 
count_s_c=0
sh=0
tot = np.zeros(len(varnames))
print(len(filenames_cycle))
for inc,ncname in enumerate(filenames_cycle) :
    i= inc % nranks
    if i == rank :
        try:
            f = nc.Dataset(ncname)
        except:
            pass
        for it,tname in enumerate(varnames):    #keep same order of xml file
            var = f.variables[tname][:,:,:,:]
            if len(var[:,0,0,0])!=lenght:
                break
            else :
                output_cycles[it,:] += var[:,0,0,0]
                var_cycles[inc,it] = variation(var[-int(lenght/5):,0,0,0])
if rank == 0:
    val = np.zeros((len(varnames),lenght))
    output_cycles_global = np.copy(output_cycles)
    val1 = np.zeros((len(filenames_cycle),len(varnames))) 
    var_cycles_global = np.copy(var_cycles)
    count_c_global = count_s_c
    val2=np.zeros(1)
    for inc in range(nranks) :
        i= inc % nranks
        print(i, flush=True)
        if i != 0:
            comm.Recv( [val2[:], MPI.DOUBLE], source=i, tag=1 )
            comm.Recv( [val[:,:], MPI.DOUBLE], source=i, tag=3 )
            output_cycles_global[:,:]+=val[:,:]
            count_c_global += val2[:]
    print("collecting shannon", flush=True)
    for ipc,ncname in enumerate(filenames_cycle):
        i = ipc % nranks
        if i!=0: 
            comm.Recv( [val1[:], MPI.DOUBLE], source=i, tag=2 )
            var_cycles_global[ipc,:] = val1[ipc,:]
else :
        val1=np.zeros(1)
        val2=np.zeros(1)
        val2[0]=float(count_s_c)
        comm.Send( [val2, MPI.DOUBLE], dest=0, tag=1 )
        comm.Send( [output_cycles[:,:], MPI.DOUBLE], dest=0, tag=3 )
        for ipc,ncname in enumerate(filenames_cycle) :
            i= ipc % nranks
            if i== rank :
                comm.Send( [var_cycles[ipc,:], MPI.DOUBLE], dest=0, tag=2 )
print('End of the loop',flush=True)
if rank == 0 :
    norm = 1.0 / len(filenames_cycle)
    mean_var_cycles = np.mean(var_cycles_global,0)
    print('collecting indexes',flush= True)
    print('pickling',flush=True)
    pkname = 'var.pickle'
    outfile = open(pkname,'wb')
    out_dict = {'C': mean_var_cycles, 'SC': var_cycles_global}#, 'I': idx}
    pickle.dump(out_dict,outfile)
    outfile.close()





