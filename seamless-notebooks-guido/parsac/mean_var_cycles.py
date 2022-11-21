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
filenames = glob.glob('/g100_scratch/userexternal/gocchipi/BFM_TOTAL_0dot01/*/*.nc')
filenames.sort(key=lambda x: int(''.join(filter(str.isdigit, x)))) #sort by number

indexes = np.zeros(len(filenames))
for iv,var in enumerate(filenames):
     indexes[iv] = int(var.rsplit('/')[-2])

#interested variables names
varnames = ["B1_c","P1_c","P2_c","P3_c","P4_c","Z5_c","Z6_c","Z3_c","Z4_c"]
#get the cycling simulation for small N1p

pknamein = "bfm_sensitivity_highcv.pickle"
pknameout = 'bfm_sensitivity_high_script.pickle'
#pknameout = 'bfm_sensitivity_total.pickle'
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
#outputs = [outputs[int(i)] for i in indexes]
#get target names

mydoc = minidom.parse('bfm_sensitivity.xml')
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)
outputs_living = np.copy(outputs).tolist()
target_living = np.copy(target).tolist()
#remove 05 03c 03h R3c O2_o
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

#only living species
right_index_living = []
varnames = ["B1_c","P1_c","P2_c","P3_c","P4_c","Z5_c","Z6_c","Z3_c","Z4_c"]
for i,name in enumerate(target_living):
    for vname in varnames:
        if name.find(vname) == 0 :
            right_index_living.append(i)
for i in range(len(outputs_living)):  #remove O5
         outputs_living[i] = [outputs_living[i][ind] for ind in sorted(right_index_living, reverse=True)]
target_living = [target_living[ind] for ind in sorted(right_index_living, reverse=True)]

count_living = np.zeros(len(inputs))
for o in range(len(outputs_living)):
    for i in range(int(len(outputs_living[o])/2)):
        if outputs_living[o][2*i]+outputs_living[o][2*i+1]==2:
            count_living[o]=1
        else :
            pass
lll = 0
for l in count_living:
    if l==1:
        lll+=1
#select directories with cycling behaviour
cyclesind =[]
for i in range(len(inputs)):
    if count[i]==1:# and inputs[i][N1p]<0.35 : 
        cyclesind.append(i)
#select directories with cycling behaviour for small N1p (<0.35)
cyclesind_low =[]
for i in range(len(inputs)):
    if count[i]==1 and inputs[i][N1p]<0.35 :
        cyclesind_low.append(i)
#corrsponding indexes between low nut and all nut
id_nut = np.zeros(len(cyclesind_low))
for i,v in enumerate(cyclesind_low):
    for j,u in enumerate(cyclesind):
        if u==v:
            id_nut[i]=j

#select directories with non-cycling behaviour for large N1p (>0.35)
linind = []
for i in range(len(inputs)):
    if count[i]==0:# and inputs[i][N1p]>0.35 :
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
            #idx = f.variables.keys().index(tname)
            var = f.variables[tname][:,:,:,:]
            #for key,var in zip(f.variables.keys(),f.variables.values()):
            #    if tname == key:
            if len(var[:,0,0,0])!=lenght:
                break
            else :
                output_cycles[it,:] += var[:,0,0,0]
                var_cycles[inc,it] = variation(var[-int(lenght/5):,0,0,0])
        #if np.isnan(np.sum(p_cycles/np.sum(p_cycles)*np.log(p_cycles/np.sum(p_cycles)))):
        #        count_s_c +=1
                #break
        #else :
        #        shannon_cycles += -np.sum(p_cycles/np.sum(p_cycles)*np.log(p_cycles/np.sum(p_cycles)))
            
            #for ik,key in enumerate(f.variables.keys()):
            #    if tname == key :
            #        break
            #for iv,var in enumerate(f.variables.values()):
            #     if iv == ik :
            #         output_cycles[it,:] += var[:,0,0,0]
            #         p_cycles[it] = np.mean(var[-int(lenght/10):,0,0,0])
            #if np.isnan(np.sum(p_cycles/np.sum(p_cycles)*np.log(p_cycles/np.sum(p_cycles)))):
            #    count_s_c +=1
            #else :
            #    shannon_cycles += -np.sum(p_cycles/np.sum(p_cycles)*np.log(p_cycles/np.sum(p_cycles)))
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
            #comm.Recv( [val1[:], MPI.DOUBLE], source=i, tag=2 )
            comm.Recv( [val[:,:], MPI.DOUBLE], source=i, tag=3 )
            output_cycles_global[:,:]+=val[:,:]
            #shannon_cycles_global += val1[:]
            count_c_global += val2[:]
    print("collecting shannon", flush=True)
    for ipc,ncname in enumerate(filenames_cycle):
        i = ipc % nranks
        if i!=0: 
            comm.Recv( [val1[:], MPI.DOUBLE], source=i, tag=2 )
            var_cycles_global[ipc,:] = val1[ipc,:]
            #print("inc: "+str(inc),flush=True)
            #print(val[:],flush=True)
else :
        val1=np.zeros(1)
        val2=np.zeros(1)
        #val1=shannon_cycles
        val2[0]=float(count_s_c)
        #comm.Send( [val1, MPI.DOUBLE], dest=0, tag=2 )
        comm.Send( [val2, MPI.DOUBLE], dest=0, tag=1 )
        comm.Send( [output_cycles[:,:], MPI.DOUBLE], dest=0, tag=3 )
        for ipc,ncname in enumerate(filenames_cycle) :
            i= ipc % nranks
            if i== rank :
                comm.Send( [var_cycles[ipc,:], MPI.DOUBLE], dest=0, tag=2 )
#       print(output_cycles[inc,:],flush=True)
print('End of the loop',flush=True)
##print('starting non cycling')
## non-cycling trajectories
#output_lin = np.zeros((len(varnames),lenght))
#filenames_lin = [filenames[i] for i in linind]
#count = 0
#count_s_nc=0
#p_lin = np.zeros(len(varnames))
#shannon_lin = np.zeros((len(filenames_lin),1)) 
#sh=0
#tot = np.zeros(len(varnames))
#for inc,ncname in enumerate(filenames_lin) :
#    i= inc % nranks
#    if i == rank :
#        f = nc.Dataset(ncname)
#        for it,tname in enumerate(varnames):    #keep same order of xml file
#            #for key,var in zip(f.variables.keys(),f.variables.values()):
#            #    if tname == key:
#            var = f.variables[tname][:,:,:,:]
#            if np.isnan(np.sum(var[:,0,0,0])):
#                  count += 1
#                  break
#            else :
#                  output_lin[it,:] += var[:,0,0,0]
#                  p_lin[it] = np.mean(var[-int(lenght/10):,0,0,0])
#                  if p_lin[it]<0:
#                      p_lin[it]=0
#        for ip,pv in enumerate(p_lin):
#            if pv ==0:
#                tot[ip]=0
#            else :
#                tot[ip]= pv/np.sum(p_lin)*np.log(pv/np.sum(p_lin))
#        sh = np.sum(tot)
#        if -sh <= np.log(9):
#            shannon_lin[inc,0] = -sh
#        else:
#            count_s_nc =+1
#
#        #if np.isnan(np.sum(p_lin/np.sum(p_lin)*np.log(p_lin/np.sum(p_lin)))):
#        #        count_s_nc +=1
#                #break
#        #else :
#        #        shannon_lin += -np.sum(p_lin/np.sum(p_lin)*np.log(p_lin/np.sum(p_lin)))
#            
#            #for ik,key in enumerate(f.variables.keys()):
#            #    if tname == key :
#            #        break
#            #for iv,var in enumerate(f.variables.values()):
#            #     if iv == ik :
#            #         if np.isnan(np.sum(var[:,0,0,0])):
#            #             count += 1
#            #         else :
#            #             output_lin[it,:] += var[:,0,0,0] 
#            #             p_lin[it] = np.mean(var[-int(lenght/10):,0,0,0])
#            #if np.isnan(np.sum(p_lin/np.sum(p_lin)*np.log(p_lin/np.sum(p_lin)))):
#            #    count_s_nc +=1
#            #else :
#            #    shannon_lin += -np.sum(p_lin/np.sum(p_lin)*np.log(p_lin/np.sum(p_lin)))
#if rank == 0:
#    val = np.zeros((len(varnames),lenght))
#    output_lin_global = np.copy(output_lin)
#    val1 = np.zeros(1) 
#    val2 = np.zeros(1)
#    val3 = np.zeros(1)
#    count_global = count
#    count_snc_global = count_s_nc
#    shannon_lin_global = shannon_lin
#    for inc in range(nranks) :
#        i= inc % nranks
#        print(i, flush=True)
#        if i != 0:
#           # comm.Recv( [val1[:], MPI.DOUBLE], source=i, tag=1 )
#            comm.Recv( [val2[:], MPI.DOUBLE], source=i, tag=2 )
#            comm.Recv( [val3[:], MPI.DOUBLE], source=i, tag=3 )
#            comm.Recv( [val[:,:], MPI.DOUBLE], source=i, tag=4 )
#            output_lin_global[:,:]+=val[:,:]
#           # shannon_lin_global += val1[:]
#            count_global += val2[:]
#            count_snc_global += val3[:]
#            #print("inc: "+str(inc),flush=True)
#            #print(val[:],flush=True)
#    print('collecting shannon',flush=True)
#    for inp,ncname in enumerate(filenames_lin):
#        i= inp % nranks
#        if i!=0:
#            comm.Recv( [val1[:], MPI.DOUBLE], source=i, tag=1 )
#            shannon_lin_global[inp,0] = val1[0]
#
#else :
#
#        val3 = np.zeros(1)
#        val2 = np.zeros(1)
#        #val1 = np.zeros(1)
#        #val1=shannon_lin
#        val2[0]=float(count)
#        val3[0]=float(count_s_nc)
#        #comm.Send( [val1, MPI.DOUBLE], dest=0, tag=1 )
#        comm.Send( [val2, MPI.DOUBLE], dest=0, tag=2 )
#        comm.Send( [val3, MPI.DOUBLE], dest=0, tag=3 )
#        comm.Send( [output_lin[:,:], MPI.DOUBLE], dest=0, tag=4 )
#        for ipc,ncname in enumerate(filenames_lin):
#            i = ipc % nranks
#            if i== rank :
#                comm.Send( [shannon_lin[ipc,0], MPI.DOUBLE], dest=0, tag=1 )
##       print(output_lin[inc,:],flush=True)
#print('End of the loop',flush=True)
#
#compute means
if rank == 0 :
    norm = 1.0 / len(filenames_cycle)
#    mean_cycles = output_cycles_global * norm
#    print(count_c_global,len(filenames_cycle))
    #norm_s_c = 1.0 / (len(filenames_cycle)-count_c_global)
    #mean_shannon_cycles = shannon_cycles_global * norm#_s_c
    mean_var_cycles = np.mean(var_cycles_global,0)
#    print(mean_shannon_cycles)#, shannon_cycles_global/len(filenames_cycle))
#    print(count_global,count_snc_global,len(filenames_lin))
#    norm1 = 1.0 / (len(filenames_lin)-count_global)
#    mean_lin = output_lin_global * norm1
##    norm_s_nc = 1.0 / (len(filenames_cycle)-count_global-count_snc_global)
#    mean_shannon_lin = np.mean(shannon_lin_global)# * norm1#_s_nc
#    print(mean_shannon_lin, shannon_lin_global/len(filenames_lin))
    print('collecting indexes',flush= True)
    idx = 0 
    for iv,v in enumerate(var_cycles_global):
        if v.any() > 0.05: #collect trajectories with CV bigger than 0.05
                idx+=1
    print('hich_cv',idx)
    print('len non-stat all',len(filenames_cycle))
    print('len non-stat living',lll)
    print('pickling',flush=True)
    pkname = 'cv.pickle'
    outfile = open(pkname,'wb')
    out_dict = {'C': mean_var_cycles, 'SC': var_cycles_global, 'I': id_nut}
    pickle.dump(out_dict,outfile)
    outfile.close()





