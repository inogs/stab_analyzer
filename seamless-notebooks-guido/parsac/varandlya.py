import pickle
import glob
import numpy as np
import netCDF4 as nc
from xml.dom import minidom
from scipy.stats import variation 
import lyapunov

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
outputs = [outputs[int(i)] for i in indexes]
#get target names

mydoc = minidom.parse('bfm_sensitivity.xml')
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)

#get parameters name from the xml
items = mydoc.getElementsByTagName('parameter')
parameters = []
for i in range(len(items)):
    string = items[i].attributes['variable'].value
    if i<200:
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
    else :
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])

#count cycles        
count = np.zeros(len(inputs))
for o in range(len(outputs)):
    for i in range(int(len(outputs[o])/2)):
        if outputs[o][2*i]+outputs[o][2*i+1]==2:
       # if outputs[o][2*i]==1:
            count[o]=1
        else :
            pass
#select directories with cycling behaviour
cyclesind =[]
for i in range(len(inputs)):
    if count[i]==1:# and inputs[i][N1p]<0.35 : 
        cyclesind.append(i)

#select directories with non-cycling behaviour for large N1p (>0.35)
linind = []
for i in range(len(inputs)):
    if count[i]==0:# and inputs[i][N1p]>0.35 :
        linind.append(i)
#lenght of trajectories
lenght = 10000 
# cycling trajectories
filenames_cycle = [filenames[i] for i in cyclesind]
l_cycles = np.zeros((len(filenames_cycle),len(varnames))) 
count_s_c=0
tot = np.zeros(len(varnames))

#restrict to not done folders
txtnames = glob.glob('/g100_work/icei_Rosati/seamless-notebooks-guido/parsac/test/*.txt')
txtnames.sort(key=lambda x: int(''.join(filter(str.isdigit, x)))) #sort by number
txtnames = [(v.rsplit('/')[-1]).rsplit('.')[-2] for v in txtnames]
txtnames = ['/g100_scratch/userexternal/gocchipi/BFM_TOTAL/'+v+'/result.nc' for v in txtnames]
newnames = [x for x in filenames_cycle if x not in txtnames]

for inc,ncname in enumerate(newnames[:1000]) :
    i= inc % nranks
    if i == rank :
        try:
            f = nc.Dataset(ncname)
        except:
            pass
        filename = 'test/'+ncname.rsplit('/')[-2]+'.txt'
        with open(filename, "w+") as o:
            o.write('')
        for it,tname in enumerate(varnames):    #keep same order of xml file
                var = f.variables[tname][:,:,:,:]
                lyap = lyapunov.LYAP(var[-int(9000):,0,0,0])
                if len(var[:,0,0,0])!=lenght:
                    break
                else :
                    try:
                        l_cycles[inc,it] = lyap.lyap_e(dt=3650./10000.)
                        with open(filename, "a") as o:
                            o.write(str(l_cycles[inc,it])+'\n')
                        print('created '+filename,flush=True)
                    except:
                        print('error: '+ncname+' '+tname,flush=True)
                        count_s_c +=1
                        pass
#comm.Barrier()
print("collecting ranks",flush=True)
if rank == 0:
    l_cycles_global = np.copy(l_cycles)
    count_c_global = count_s_c
    val2=np.zeros(1)
    for inc in range(nranks) :
        i= inc % nranks
        print(i, flush=True)
        if i != 0:
            comm.Recv( [val2[:], MPI.DOUBLE], source=i, tag=1 )
            count_c_global += val2[:]
    for ipc,ncname in enumerate(filenames_cycle):
        i = ipc % nranks
        if i!=0: 
            lya = comm.recv(source=i, tag=4 )
            l_cycles_global[int(lya['idx']),:] = lya['data']
else :
        val2=np.zeros(1)
        val2[0]=float(count_s_c)
        comm.Send( [val2, MPI.DOUBLE], dest=0, tag=1 )
        for ipc,ncname in enumerate(filenames_cycle) :
            i= ipc % nranks
            lya = {'idx':ipc, 'data': l_cycles[ipc,:,]}
            if i== rank :
                comm.send( lya, dest=0, tag=4 )
comm.Barrier()
#       print(output_cycles[inc,:],flush=True)
print('End of the loop',flush=True)
#compute means
if rank == 0 :
    print('pickling',flush=True)
    pkname = 'var.pickle'
    outfile = open(pkname,'wb')
    out_dict = {'L': l_cycles_global}
    pickle.dump(out_dict,outfile)
    outfile.close()
    print('errors: ',count_c_global)
    print('not-steady: ', len(filenames_cycle))




