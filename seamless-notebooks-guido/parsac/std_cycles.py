import pickle
import glob
import numpy as np
import netCDF4 as nc
from xml.dom import minidom
import sys
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
key = 'O5'
keyOc= 'O3_c'
keyOh= 'O3h_h'
keyR3c = 'R3_c'
keyO2 = 'O2_o'
false_index = []
for i,name in enumerate(target):
    if name.find(keyO2) !=-1 or name.find(key) !=-1 or name.find(keyOc) !=-1 or name.find(keyOh) !=-1 or name.find(keyR3c) !=-1:
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

## cycling trajectories
#output_cycles = np.zeros((len(varnames),lenght))
#filenames_cycle = [filenames[i] for i in cyclesind]
#print(len(filenames_cycle))
#for inc,ncname in enumerate(filenames_cycle) :
#    i= inc % nranks
#    if i == rank :
#        f = nc.Dataset(ncname)
#        for it,tname in enumerate(varnames):    #keep same order of xml file
#            for ik,key in enumerate(f.variables.keys()):
#                if tname == key :
#                    break
#            for iv,var in enumerate(f.variables.values()):
#                 if iv == ik :
#                     output_cycles[it,:] += var[:,0,0,0]
#if rank == 0:
#    val = np.zeros((len(varnames),lenght))
#    output_cycles_global = np.copy(output_cycles)
#    for inc in range(nranks) :
#        i= inc % nranks
#        print(i, flush=True)
#        if i != 0:
##           comm.Recv( [idx_inc, MPI.INT], source=i, tag=1 )
#            #comm.Recv( [idx_it, MPI.INT], source=i, tag=2 )
#            comm.Recv( [val[:,:], MPI.DOUBLE], source=i, tag=3 )
#            output_cycles_global[:,:]+=val[:,:]
#            #print("inc: "+str(inc),flush=True)
#            #print(val[:],flush=True)
#else :
##       comm.Send( [inc, MPI.INT], dest=0, tag=1 )
#        comm.Send( [output_cycles[:,:], MPI.DOUBLE], dest=0, tag=3 )
##       print(output_cycles[inc,:],flush=True)
#print('End of the loop',flush=True)
#
#
## non-cycling trajectories
#output_lin = np.zeros((len(varnames),lenght))
#filenames_lin = [filenames[i] for i in linind]
#count = 0
#for inc,ncname in enumerate(filenames_lin) :
#    i= inc % nranks
#    if i == rank :
#        f = nc.Dataset(ncname)
#        for it,tname in enumerate(varnames):    #keep same order of xml file
#            for ik,key in enumerate(f.variables.keys()):
#                if tname == key :
#                    break
#            for iv,var in enumerate(f.variables.values()):
#                 if iv == ik :
#                     if np.isnan(np.sum(var[:,0,0,0])):
#                         count += 1
#                     else :
#                         output_lin[it,:] += var[:,0,0,0] 
#if rank == 0:
#    val = np.zeros((len(varnames),lenght))
#    output_lin_global = np.copy(output_lin)
#    for inc in range(nranks) :
#        i= inc % nranks
#        print(i, flush=True)
#        if i != 0:
##           comm.Recv( [idx_inc, MPI.INT], source=i, tag=1 )
#            #comm.Recv( [idx_it, MPI.INT], source=i, tag=2 )
#            comm.Recv( [val[:,:], MPI.DOUBLE], source=i, tag=3 )
#            output_lin_global[:,:]+=val[:,:]
#            #print("inc: "+str(inc),flush=True)
#            #print(val[:],flush=True)
#else :
##       comm.Send( [inc, MPI.INT], dest=0, tag=1 )
#        comm.Send( [output_lin[:,:], MPI.DOUBLE], dest=0, tag=3 )
##       print(output_lin[inc,:],flush=True)
#print('End of the loop',flush=True)
#
##compute means
#if rank == 0 :
#    norm = 1.0 / len(filenames_cycle)
#    mean_cycles = output_cycles_global * norm
#    print(count,len(filenames_lin))
#    norm1 = 1.0 / (len(filenames_lin)-count)
#    mean_lin = output_lin_global * norm1
#print('inving means between ranks', flush=True)
#if rank == 0:
#    val1=np.copy(mean_cycles)
#    val2=np.copy(mean_lin)
#    comm.Send( [val1, MPI.DOUBLE], dest=0, tag=2 )
#    comm.Send( [val2, MPI.DOUBLE], dest=0, tag=1 )
#    for inc in range(nranks) :
#        i= inc % nranks
#        print(i, flush=True)
#        if i != 0:
#            comm.Send( [val1, MPI.DOUBLE], dest=i, tag=1 )
#            comm.Send( [val2, MPI.DOUBLE], dest=i, tag=2 )
#else :
#        mean_cycles = np.zeros(lenght)
#        mean_lin = np.zeros(lenght)
#        val1=np.zeros(lenght)
#        val2=np.zeros(lenght)
#        comm.Recv( [val1[:], MPI.DOUBLE], source=0, tag=1 )
#        comm.Recv( [val2[:], MPI.DOUBLE], source=0, tag=2 )
#        for inc in range(nranks) :
#            i= inc % nranks
#            print(i, flush=True)
#            if i != 0:
#                mean_cycles = val1[:]
#                mean_lin = val2[:]
#
#print('start cycles', flush=True)

pkname = 'means.pickle'
infile = open(pkname,'rb')
new_dict = pickle.load(infile)
infile.close()
mean_cycles = new_dict.get('C')
mean_lin = new_dict.get('NC')
# cycling trajectories
output_cycles = np.zeros((len(varnames),lenght))
filenames_cycle = [filenames[i] for i in cyclesind]
print(len(filenames_cycle))
for inc,ncname in enumerate(filenames_cycle) :
    i= inc % nranks
    if i == rank :
        f = nc.Dataset(ncname)
        for it,tname in enumerate(varnames):    #keep same order of xml file
            var = f.variables[tname][:,:,:,:]         
            if len(var[:,0,0,0]) != lenght:
                break
            else :
                output_cycles[it,:] += np.square(var[:,0,0,0]-mean_cycles[it])
if rank == 0:
    val = np.zeros((len(varnames),lenght))
    output_cycles_global = np.copy(output_cycles)
    for inc in range(nranks) :
        i= inc % nranks
        print(i, flush=True)
        if i != 0:
#           comm.Recv( [idx_inc, MPI.INT], source=i, tag=1 )
            #comm.Recv( [idx_it, MPI.INT], source=i, tag=2 )
            comm.Recv( [val[:,:], MPI.DOUBLE], source=i, tag=3 )
            output_cycles_global[:,:]+=val[:,:]
            #print("inc: "+str(inc),flush=True)
            #print(val[:],flush=True)
else :
#       comm.Send( [inc, MPI.INT], dest=0, tag=1 )
        comm.Send( [output_cycles[:,:], MPI.DOUBLE], dest=0, tag=3 )
#       print(output_cycles[inc,:],flush=True)
print('End of the loop',flush=True)

print('start non-cycles', flush=True)
# non-cycling trajectories
output_lin = np.zeros((len(varnames),lenght))
filenames_lin = [filenames[i] for i in linind]
count = 0
for inc,ncname in enumerate(filenames_lin) :
    i= inc % nranks
    if i == rank :
        f = nc.Dataset(ncname)
        for it,tname in enumerate(varnames):    #keep same order of xml file
            var = f.variables[tname][:,:,:,:]
            if np.isnan(np.sum(var[:,0,0,0])):
                  count += 1
                  break
            elif len(var[:,0,0,0])!=lenght:
                break
            else :
                output_lin[it,:] += np.square(var[:,0,0,0]-mean_lin[it])
if rank == 0:
    val = np.zeros((len(varnames),lenght))
    val1 = np.zeros(1)
    count_global = count
    output_lin_global = np.copy(output_lin)
    for inc in range(nranks) :
        i= inc % nranks
        print(i, flush=True)
        if i != 0:
#           comm.Recv( [idx_inc, MPI.INT], source=i, tag=1 )
            #comm.Recv( [idx_it, MPI.INT], source=i, tag=2 )
            comm.Recv( [val[:,:], MPI.DOUBLE], source=i, tag=3 )
            comm.Recv( [val1[:], MPI.DOUBLE], source=i, tag=2 )
            output_lin_global[:,:]+=val[:,:]
            count_global += val1[:]
            #print("inc: "+str(inc),flush=True)
            #print(val[:],flush=True)
else :
        val1 = np.zeros(1)
        val1[0] = float(count)
#       comm.Send( [inc, MPI.INT], dest=0, tag=1 )
        comm.Send( [val1, MPI.DOUBLE], dest=0, tag=2 )
        comm.Send( [output_lin[:,:], MPI.DOUBLE], dest=0, tag=3 )
#       print(output_lin[inc,:],flush=True)
print('End of the loop',flush=True)
print('Writing pickle', flush=True)
#compute means
if rank == 0 :
    norm = 1.0 / len(filenames_cycle)
    std_cycles = np.sqrt(output_cycles_global * norm)
    print(len(filenames_cycle))
    norm1 = 1.0 / (len(filenames_lin)-count_global)
    std_lin = np.sqrt(output_lin_global * norm1)
    pkname = 'std.pickle'
    outfile = open(pkname,'wb')
    out_dict = {'C': std_cycles, 'NC' : std_lin}
    pickle.dump(out_dict,outfile)
    outfile.close()
