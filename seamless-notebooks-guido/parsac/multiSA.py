import pickle
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap



#function to give the rank to the sensitivities
def changeArr(input1):

    # Copy input array into newArray
    newArray = input1.copy()

    # Sort newArray[] in ascending order
    newArray.sort()

    # Dictionary to store the rank of
    # the array element
    ranks = {}

    rank = 1

    for index in range(len(newArray)):
        element = newArray[index];

        # Update rank of element
        if element not in ranks:
            ranks[element] = rank
            rank += 1

    # Assign ranks to elements
    for index in range(len(input1)):
        element = input1[index]
        input1[index] = ranks[input1[index]]





filename1='bfm_sensitivity_fast.analyze.pickle'   #fast
filename2='rma_sensitivity.analyze.pickle'
filename3='twotwo_sensitivity.analyze.pickle'   #sobol
filename4='bfm_sensitivity_fast_mvr.analyze.pickle'  #must be mvr to get correct names of parameters and outputs

infile = open(filename1,'rb')
new_dict1 = pickle.load(infile)
infile.close()

infile = open(filename2,'rb')
new_dict2 = pickle.load(infile)
infile.close()

infile = open(filename3,'rb')
new_dict3 = pickle.load(infile)
infile.close()

infile = open(filename4,'rb')
new_dict4 = pickle.load(infile)
infile.close()


substringV = 'var'
substringF = 'peak'
Noutputs = int(len(new_dict1.keys())/2)
#var = np.zeros(Noutputs)
var =[]
four = []

#file1

for key in new_dict1.keys():
    if key.find(substringV) != -1 :
        a = new_dict1.get(key)
        for name in a.keys() :
            if name == 'S1' or name == 'beta':
                var.append(np.abs(a.get(name)).tolist())
    elif key.find(substringF) != -1 :
        a = new_dict1.get(key)
        for name in a.keys() :
            if name == 'S1' or name == 'beta':
                four.append(np.abs(a.get(name)).tolist())
SA1 = var.copy()
#give the rank to the sensitivities
for i in range(len(var)):
    changeArr(var[i])   
    changeArr(four[i])
    SA1[i] = [x+y for x,y in zip(var[i],four[i])]

#file2
#var.clear()
#four.clear()
#
#for key in new_dict2.keys():
#    if key.find(substringV) != -1 :
#        a = new_dict2.get(key)
#        for name in a.keys() :
#            if name == 'S1' or name == 'beta':
#                var.append(np.abs(a.get(name)).tolist())
#    elif key.find(substringF) != -1 :
#        a = new_dict2.get(key)
#        for name in a.keys() :
#            if name == 'S1' or name == 'beta':
#                four.append(np.abs(a.get(name)).tolist())
#SA2 = var.copy()
#give the rank to the sensitivities
#for i in range(len(var)):
#    changeArr(var[i])
#    changeArr(four[i])
#    SA2[i] = [x+y for x,y in zip(var[i],four[i])]

#file3  sobol
#var.clear()
#four.clear()
#
#for key in new_dict3.keys():
#    if key.find(substringV) != -1 :
#        a = new_dict3.get(key)
#        for name in a.keys() :
#            if name == 'ST' or name == 'beta':
#                var.append(np.abs(a.get(name)).tolist())
#    elif key.find(substringF) != -1 :
#        a = new_dict3.get(key)
#        for name in a.keys() :
#            if name == 'ST' or name == 'beta':
#                four.append(np.abs(a.get(name)).tolist())
#SA3 = var.copy()
#give the rank to the sensitivities
#for i in range(len(var)):
#    changeArr(var[i])
#    changeArr(four[i])
#    SA3[i] = [x+y for x,y in zip(var[i],four[i])]

#file4
var.clear()
four.clear()

for key in new_dict4.keys():
    if key.find(substringV) != -1 :
        a = new_dict4.get(key)
        for name in a.keys() :
            if name == 'ST' or name == 'beta':
                var.append(np.abs(a.get(name)).tolist())
    elif key.find(substringF) != -1 :
        a = new_dict4.get(key)
        for name in a.keys() :
            if name == 'ST' or name == 'beta':
                four.append(np.abs(a.get(name)).tolist())
SA4 = var.copy()
#give the rank to the sensitivities
for i in range(len(var)):
    changeArr(var[i])
    changeArr(four[i])
    SA4[i] = [x+y for x,y in zip(var[i],four[i])]


#sum the sensitivities from the different methods
SA = SA1.copy()
for i in range(len(SA1)):
    SA[i] = [x+y for x,y in zip(SA1[i],SA4[i])]#,SA3[i],SA2[i])] mancano z,w
    changeArr(SA[i])

#define the parameters names file4 MUST be MVR
parameters = []#np.zeros(len(a.get('names')), dtype=str)
for i,string in enumerate(a.get('names')):
    parameters.append(string.rsplit('/')[-1])

#define the output names from the file result.nc

fn = '/g100_work/OGS21_PRACE_P/seamless-notebooks-guido/setups/ogs/result.nc'
ds = nc.Dataset(fn)
var = ds.variables
var.pop('z')
var.pop('lat')
var.pop('lon')
var.pop('time')
output = []
for key in var.keys():
    if key.find(substringV) != -1 :
        break
    else :
        output.append(key)

#colorplot

#Nrank = len(SA[0])
#x = np.arange(1, Nrank+1, 1)
#y = np.arange(1, len(output)+1, 1)
#plt.hist2d(x, y, bins=Nrank, weights=SA)
#plt.xticks(x, parameters) #parameters' names on x axis
#plt.yticks(y, output) #outputs' names on y axis
#plt.show()

#Sum the sensitivities of all parameters to obtain the most relevant ones for the whole model
SA_tot = SA[0].copy()
for i in range(len(SA)):
    SA_tot = [x+y for x,y in zip(SA_tot,SA[i])]

#Nparam more sensible parameters
Nparam=30
SA_tot = np.array(SA_tot)
SA = np.array(SA)
indecies = np.zeros(len(parameters), dtype=int) 
param = [] 
param_lim = []
SA_lim = []
indecies = np.argsort(SA_tot)
param = [parameters[l] for l in indecies]
for i in range(len(SA)):
    SA[i] = SA[i][indecies]
for i in range(len(SA)):
    SA_lim.append(SA[i][-Nparam:])
param_lim=param[-Nparam:]
#colorplot
Nrank = len(SA_lim[0])
fig, axs = plt.subplots(len(SA_lim),1,figsize=(13,50))
fig.tight_layout()
x = np.arange(1, Nrank+1, 1)
y = np.ones(len(x))
for isa,sa in enumerate(SA_lim) :
    im = axs[isa].scatter(x, y, c=sa)
    axs[isa].set_xticks(x, minor=False) #parameters' names on x axis
    axs[isa].set_xticklabels(param_lim, fontdict=None, minor=False)
    axs[isa].set_ylabel(output[isa])
for ax in fig.get_axes():
    ax.tick_params(axis='both', which='major', labelsize=5)
    ax.tick_params(axis='both', which='minor', labelsize=5)

cb_ax = fig.add_axes([0.97, 0.10, 0.002, 0.8])
cbar = fig.colorbar(im, cax=cb_ax)
#set the colorbar ticks and tick labels
cbar.set_ticks(np.arange(0, 1.1, 0.5))
cbar.set_ticklabels(['low', 'medium', 'high'])

fig.savefig('SA30.png', format='png',dpi=150)



