import pickle
import glob
import numpy as np
import netCDF4 as nc
from xml.dom import minidom

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
filenames = glob.glob('/g100_scratch/userexternal/gocchipi/FUSSMAN_ZOOM/*/*.nc')
filenames.sort(key=lambda x: int(''.join(filter(str.isdigit, x)))) #sort by number

indexes = np.zeros(len(filenames))
for iv,var in enumerate(filenames):
     indexes[iv] = int(var.rsplit('/')[-2])

#get target names

mydoc = minidom.parse('fuss_sensitivity.xml')
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)

#get the target values:
print('nranks= '+str(nranks), flush=True)
print('myrank= '+str(rank), flush=True)
print(len(filenames), flush=True)
output = np.zeros((len(filenames),len(target)))
#for inc,ncname in enumerate(filenames[rank::nranks]) :
for inc,ncname in enumerate(filenames) :
    i= inc % nranks
    if i == rank :
        try:
            f = nc.Dataset(ncname)
        except:
            pass   
        for it,tname in enumerate(target):    #keep same order of xml file
            for ik,key in enumerate(f.variables.keys()):
                if tname == key :
                    break
            for iv,var in enumerate(f.variables.values()):
                 if iv == ik :
                     output[inc,it] = var[:]
if rank == 0:
    val = np.zeros(len(target))
    for inc,ncname in enumerate(filenames) :
        i= inc % nranks
        print(i, flush=True)
        if i != 0:
            out = comm.recv( source=i, tag=1 )
            output[int(out['idx']),:]=out['data']
else :
    for inc,ncname in enumerate(filenames) :
        i = inc % nranks
        out = {'idx':inc, 'data':output[inc,:]}
        if i == rank :
            comm.send( out, dest=0, tag=1 )
print('End of the loop',flush=True)
#add to the pickle file
if rank ==0:
    print('Printing pickle...',flush=True)
    pkname = 'fuss_sensitivity_script.pickle'
    #infile = open(pkname,'rb')
    #new_dict = pickle.load(infile)
    #infile.close()
    #inputs = new_dict.get('X')
    outfile = open(pkname,'wb')
    out_dict = {'I':indexes, 'Y' : output}
    pickle.dump(out_dict,outfile)
    outfile.close()





