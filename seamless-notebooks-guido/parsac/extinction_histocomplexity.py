import pickle
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import argparse
from xml.dom import minidom
import glob

dashList = [(0,(1, 10)),(0,(1, 1)),(0,(5,1)),(0,(3,5,1,5)),(0,(3,1,1,1,1,1))] 
color = ["black","maroon","cornflowerblue","olivedrab","turquoise"]

####FUSMANN 5D
dd=0
#get parameters name from the xml
mydoc = minidom.parse('bfm_sensitivity_5D.xml')
items = mydoc.getElementsByTagName('parameter')
parameters = []
for i in range(len(items)):
    string = items[i].attributes['variable'].value
    if i<200:
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
    else :
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
#get targets names
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)
#remove 05 03c 03h R3c

#### REMOVE ALSO O2O

key = 'O5'
keyOc= 'O3_c'
keyO2o= 'O2_o'
keyOh= 'O3h_h'
keyR3c = 'R3_c'
false_index = []
for i,name in enumerate(target):
    if name.find(key) == 0 or name.find(keyOc) == 0 or name.find(keyOh) == 0 or name.find(keyR3c) == 0 or name.find(keyO2o) == 0:
        false_index.append(i)
#target.pop()

infile = open('bfm_sensitivity_5D.pickle','rb')
new_dict = pickle.load(infile)
infile.close()
infile1 = open('bfm_sensitivity_long_extinction.pickle','rb')
new_dict1 = pickle.load(infile1)
infile1.close()
inputs = new_dict.get('X')
outputs = new_dict1.get('Y')
outputs = outputs.tolist()
#limit inputs to the folders generated
indexes = new_dict1.get('I')
kk=0
for v in indexes:
    if v!=0:
        kk+=1
print('living: ',kk)
inputs = [inputs[int(i)] for i in indexes]
#outputs = [outputs[int(i)] for i in indexes]
outputs_living = np.copy(outputs).tolist()
target_living = np.copy(target).tolist()
print('len',len(inputs))

for i in range(len(outputs)):  #remove O5
    for ind in sorted(false_index, reverse=True):
         outputs[i].pop(ind)
for ind in sorted(false_index, reverse=True):
    target.pop(ind)

count = np.zeros(len(inputs))
for o in range(len(outputs)):
    for i in range(int(len(outputs[o])/2)):
        if outputs[o][2*i]+outputs[o][2*i+1]==2:
    #if np.sum(outputs[o][1::2]) >1:
            count[o]=1
        else :
            pass#count[o]=0
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
print('non-stationary',lll)
print('long chain fraction',lll/len(outputs_living))

#histo only for living counts

count = count_living
outputs = outputs_living
target = target_living

x = np.zeros(len(inputs))
####
ll = 0
for v in count:
    if v ==1:
        ll+=1
#print('long chain non-stat rate: ', ll/len(count))
###
ind_n = []
for iv,var in enumerate(parameters):
    if var == 'N3n' or var=='N4n':
        ind_n.append(iv)

x_n=np.zeros(len(inputs))
for i in ind_n:
    for j in range(len(inputs)):
        x_n[j]+=inputs[j][i]
hist,bin_edge = np.histogram(x_n, weights=count)
hist_norm,bin_ed = np.histogram(x_n)
bins=np.zeros(len(hist))
for b in range(len(bin_edge)-1):
        bins[b]=(bin_edge[b+1]+bin_edge[b])/2

fig3, axs3 = plt.subplots(1,2)
plt.rcParams.update({'font.size': 13})
axs3 = axs3.ravel()
axs31 = axs3[1].plot(bins, hist/hist_norm,c=color[dd],label='Long chain')
#axs31 = axs3[1].plot(bins, hist/hist_norm,c='black',linestyle=dashList[dd],label='Long chain')
for i in range(len(inputs[0])):
    if parameters[i] == 'N1p':
        for j in range(len(inputs)):
            x[j]=inputs[j][i]
        hist,bin_edge = np.histogram(x, weights=count)
        hist_norm,bin_ed = np.histogram(x)
        bins=np.zeros(len(hist))
        for b in range(len(bin_edge)-1):
            bins[b]=(bin_edge[b+1]+bin_edge[b])/2
        axs30=axs3[0].plot(bins, hist/hist_norm,c=color[dd])
#       axs30=axs3[0].plot(bins, hist/hist_norm,c='black',linestyle=dashList[dd])
#axs3.set_ylim(0, 0.5)
string_name = parameters[ind_n[0]]+' + '+parameters[ind_n[1]]

### OMN. CHAIN
dd+=1
#get parameters name from the xml
mydoc = minidom.parse('bfm_sensitivity_omnchain.xml')
items = mydoc.getElementsByTagName('parameter')
parameters = []
for i in range(len(items)):
    string = items[i].attributes['variable'].value
    if i<200:
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
    else :
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
#get targets names
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)
#remove 05 03c 03h R3c

#### REMOVE ALSO O2O
#### REMOVE ALSO O2O

key = 'O5'
keyOc= 'O3_c'
keyO2o= 'O2_o'
keyOh= 'O3h_h'
keyR3c = 'R3_c'
false_index = []
for i,name in enumerate(target):
    if name.find(key) == 0 or name.find(keyOc) == 0 or name.find(keyOh) == 0 or name.find(keyR3c) == 0 or name.find(keyO2o) == 0:
        false_index.append(i)
#target.pop()

infile = open('bfm_sensitivity_omnchain_extinction.pickle','rb')
new_dict = pickle.load(infile)
infile.close()
infile1 = open('bfm_sensitivity_omnchain.pickle','rb')
new_dict1 = pickle.load(infile1)
infile1.close()
inputs = new_dict1.get('X')
outputs = new_dict.get('Y')
outputs = outputs.tolist()
#limit inputs to the folders generated
indexes = new_dict.get('I')
kk=0
for v in indexes:
    if v!=0:
        kk+=1
print('living: ',kk)
inputs = [inputs[int(i)] for i in indexes]
#outputs = [outputs[int(i)] for i in indexes]
outputs_living = np.copy(outputs).tolist()
target_living = np.copy(target).tolist()
for i in range(len(outputs)):  #remove O5
    for ind in sorted(false_index, reverse=True):
         outputs[i].pop(ind)
for ind in sorted(false_index, reverse=True):
    target.pop(ind)

count = np.zeros(len(inputs))
for o in range(len(outputs)):
    for i in range(int(len(outputs[o])/2)):
        if outputs[o][2*i]+outputs[o][2*i+1]==2:
    #if np.sum(outputs[o][1::2]) >1:
            count[o]=1
        else :
            pass#count[o]=0

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
print('omn chain non-stationary',lll)
print('fraction',lll/len(outputs_living))

#histo only for living counts

count = count_living
outputs = outputs_living
target = target_living

x = np.zeros(len(inputs))
####
ll = 0
for v in count:
    if v ==1:
        ll+=1
print('omn chain non-stat rate: ', ll/len(count))
###
ind_n = []
for iv,var in enumerate(parameters):
    if var == 'N3n' or var=='N4n':
        ind_n.append(iv)
x_n=np.zeros(len(inputs))
for i in ind_n:
    for j in range(len(inputs)):
        x_n[j]+=inputs[j][i]
hist,bin_edge = np.histogram(x_n, weights=count)
hist_norm,bin_ed = np.histogram(x_n)
bins=np.zeros(len(hist))
for b in range(len(bin_edge)-1):
        bins[b]=(bin_edge[b+1]+bin_edge[b])/2

axs31 = axs3[1].plot(bins, hist/hist_norm,c=color[dd],label='Omn. chain')
for i in range(len(inputs[0])):
    if parameters[i] == 'N1p':
        for j in range(len(inputs)):
            x[j]=inputs[j][i]
        hist,bin_edge = np.histogram(x, weights=count)
        hist_norm,bin_ed = np.histogram(x)
        bins=np.zeros(len(hist))
        for b in range(len(bin_edge)-1):
            bins[b]=(bin_edge[b+1]+bin_edge[b])/2
        axs30=axs3[0].plot(bins, hist/hist_norm,c=color[dd])
#axs3.set_ylim(0, 0.5)


####   OMNIVORY 
dd+=1
#get parameters name from the xml
mydoc = minidom.parse('bfm_sensitivity_omnivory.xml')
items = mydoc.getElementsByTagName('parameter')
parameters = []
for i in range(len(items)):
    string = items[i].attributes['variable'].value
    if i<200:
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
    else :
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
#get targets names
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)
#remove 05 03c 03h R3c

#### REMOVE ALSO O2O

key = 'O5'
keyOc= 'O3_c'
keyO2o= 'O2_o'
keyOh= 'O3h_h'
keyR3c = 'R3_c'
false_index = []
for i,name in enumerate(target):
    if name.find(key) == 0 or name.find(keyOc) == 0 or name.find(keyOh) == 0 or name.find(keyR3c) == 0 or name.find(keyO2o) == 0:
        false_index.append(i)
#target.pop()

infile = open('bfm_sensitivity_omnivory_extinction.pickle','rb')
new_dict = pickle.load(infile)
infile.close()
infile1 = open('bfm_sensitivity_omnivory.pickle','rb')
new_dict1 = pickle.load(infile1)
infile1.close()
inputs = new_dict1.get('X')
outputs = new_dict.get('Y')
outputs = outputs.tolist()
#limit inputs to the folders generated
indexes = new_dict.get('I')
kk=0
for v in indexes:
    if v!=0:
        kk+=1
print('living: ',kk)
inputs = [inputs[int(i)] for i in indexes]
outputs_living = np.copy(outputs).tolist()
target_living = np.copy(target).tolist()
for i in range(len(outputs)):  #remove O5
    for ind in sorted(false_index, reverse=True):
         outputs[i].pop(ind)
for ind in sorted(false_index, reverse=True):
    target.pop(ind)

count = np.zeros(len(inputs))
for o in range(len(outputs)):
    for i in range(int(len(outputs[o])/2)):
        if outputs[o][2*i]+outputs[o][2*i+1]==2:
    #if np.sum(outputs[o][1::2]) >1:
            count[o]=1
        else :
            pass#count[o]=0
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
print('non-stationary',lll)
print('omnivory fraction',lll/len(outputs_living))

#histo only for living counts

count = count_living
outputs = outputs_living
target = target_living

x = np.zeros(len(inputs))

ind_n = []
for iv,var in enumerate(parameters):
    if var == 'N3n' or var=='N4n':
        ind_n.append(iv)

x_n=np.zeros(len(inputs))
for i in ind_n:
    for j in range(len(inputs)):
        x_n[j]+=inputs[j][i]
hist,bin_edge = np.histogram(x_n, weights=count)
hist_norm,bin_ed = np.histogram(x_n)
bins=np.zeros(len(hist))
for b in range(len(bin_edge)-1):
        bins[b]=(bin_edge[b+1]+bin_edge[b])/2

axs31 = axs3[1].plot(bins, hist/hist_norm,c=color[dd], label='Omnivory')
for i in range(len(inputs[0])):
    if parameters[i] == 'N1p':
        for j in range(len(inputs)):
            x[j]=inputs[j][i]
        hist,bin_edge = np.histogram(x, weights=count)
        hist_norm,bin_ed = np.histogram(x)
        bins=np.zeros(len(hist))
        for b in range(len(bin_edge)-1):
            bins[b]=(bin_edge[b+1]+bin_edge[b])/2
        axs30=axs3[0].plot(bins, hist/hist_norm,c=color[dd])



### LOW GRAVITY
dd+=1
#get parameters name from the xml
mydoc = minidom.parse('bfm_sensitivity_lowgravity.xml')
items = mydoc.getElementsByTagName('parameter')
parameters = []
for i in range(len(items)):
    string = items[i].attributes['variable'].value
    if i<200:
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
    else :
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
#get targets names
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)
#remove 05 03c 03h R3c
#### REMOVE ALSO O2O

key = 'O5'
keyOc= 'O3_c'
keyO2o= 'O2_o'
keyOh= 'O3h_h'
keyR3c = 'R3_c'
false_index = []
for i,name in enumerate(target):
    if name.find(key) == 0 or name.find(keyOc) == 0 or name.find(keyOh) == 0 or name.find(keyR3c) == 0 or name.find(keyO2o) == 0:
        false_index.append(i)
#target.pop()

infile = open('bfm_sensitivity_lowgravity_extinction.pickle','rb')
new_dict = pickle.load(infile)
infile.close()
infile1 = open('bfm_sensitivity_lowgravity.pickle','rb')
new_dict1 = pickle.load(infile1)
infile1.close()

inputs = new_dict1.get('X')
outputs = new_dict.get('Y')
outputs = outputs.tolist()
#limit inputs to the folders generated
indexes = new_dict.get('I')
kk=0
for v in indexes:
    if v!=0:
        kk+=1
print('living :',kk)
inputs = [inputs[int(i)] for i in indexes]
outputs_living = np.copy(outputs).tolist()
target_living = np.copy(target).tolist()
for i in range(len(outputs)):  #remove O5
    for ind in sorted(false_index, reverse=True):
         outputs[i].pop(ind)
for ind in sorted(false_index, reverse=True):
    target.pop(ind)

count = np.zeros(len(inputs))
for o in range(len(outputs)):
    for i in range(int(len(outputs[o])/2)):
        if outputs[o][2*i]+outputs[o][2*i+1]==2:
    #if np.sum(outputs[o][1::2]) >1:
            count[o]=1
        else :
            pass#count[o]=0

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
print('non-stationary',lll)
print('low gravity fraction',lll/len(outputs_living))

#histo only for living counts

count = count_living
outputs = outputs_living
target = target_living

x = np.zeros(len(inputs))
ind_n = []
for iv,var in enumerate(parameters):
    if var == 'N3n' or var=='N4n':
        ind_n.append(iv)

x_n=np.zeros(len(inputs))
for i in ind_n:
    for j in range(len(inputs)):
        x_n[j]+=inputs[j][i]
hist,bin_edge = np.histogram(x_n, weights=count)
hist_norm,bin_ed = np.histogram(x_n)
bins=np.zeros(len(hist))
for b in range(len(bin_edge)-1):
        bins[b]=(bin_edge[b+1]+bin_edge[b])/2

axs31 = axs3[1].plot(bins, hist/hist_norm,c=color[dd],label='Low gravity')
for i in range(len(inputs[0])):
    if parameters[i] == 'N1p':
        for j in range(len(inputs)):
            x[j]=inputs[j][i]
        hist,bin_edge = np.histogram(x, weights=count)
        hist_norm,bin_ed = np.histogram(x)
        bins=np.zeros(len(hist))
        for b in range(len(bin_edge)-1):
             bins[b]=(bin_edge[b+1]+bin_edge[b])/2
        axs30=axs3[0].plot(bins, hist/hist_norm,c=color[dd])


##### BFM
dd+=1
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
#get targets names
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)
#remove 05 03c 03h R3c
#### REMOVE ALSO O2O

key = 'O5'
keyOc= 'O3_c'
keyO2o= 'O2_o'
keyOh= 'O3h_h'
keyR3c = 'R3_c'
false_index = []
for i,name in enumerate(target):
    if name.find(key) !=-1 or name.find(keyOc) !=-1 or name.find(keyOh) !=-1 or name.find(keyR3c) !=-1 or name.find(keyO2o) !=-1 :    
                false_index.append(i)
#target.pop()

infile = open('bfm_sensitivity_total_extinction.pickle','rb')
new_dict = pickle.load(infile)
infile.close()
infile1 = open('bfm_sensitivity_total_high.pickle','rb')
new_dict1 = pickle.load(infile1)
infile1.close()

inputs = new_dict1.get('X')
outputs = new_dict.get('Y')
outputs = outputs.tolist()
#limit inputs to the folders generated
indexes = new_dict.get('I')
kk=0
for v in indexes:
    if v!=0:
        kk+=1
print('living: ',kk)
inputs = [inputs[int(i)] for i in indexes]
#outputs = [outputs[int(i)] for i in indexes]
outputs_living = np.copy(outputs).tolist()
target_living = np.copy(target).tolist()
for i in range(len(outputs)):  #remove O5
    for ind in sorted(false_index, reverse=True):
         outputs[i].pop(ind)
for ind in sorted(false_index, reverse=True):
    target.pop(ind)
count = np.zeros(len(inputs))
for o in range(len(outputs)):
    for i in range(int(len(outputs[o])/2)):
        if outputs[o][2*i]+outputs[o][2*i+1]==2:
    #if np.sum(outputs[o][1::2]) >1:
            count[o]=1
        else :
            pass#count[o]=0

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
print('non-stationary',lll)
print('bfm fraction',lll/len(outputs_living))

#histo only for living counts

count = count_living
outputs = outputs_living
target = target_living

x = np.zeros(len(inputs))
ind_n = []
for iv,var in enumerate(parameters):
    if var == 'N3n' or var=='N4n':
        ind_n.append(iv)

x_n=np.zeros(len(inputs))
for i in ind_n:
    for j in range(len(inputs)):
        x_n[j]+=inputs[j][i]
hist,bin_edge = np.histogram(x_n, weights=count)
hist_norm,bin_ed = np.histogram(x_n)
bins=np.zeros(len(hist))
for b in range(len(bin_edge)-1):
        bins[b]=(bin_edge[b+1]+bin_edge[b])/2
axs31 = axs3[1].plot(bins, hist/hist_norm,c=color[dd],label='BFM')
for i in range(len(inputs[0])):
    if parameters[i] == 'N1p':
        for j in range(len(inputs)):
            x[j]=inputs[j][i]
        hist,bin_edge = np.histogram(x, weights=count)
        hist_norm,bin_ed = np.histogram(x)
        bins=np.zeros(len(hist))
        for b in range(len(bin_edge)-1):
             bins[b]=(bin_edge[b+1]+bin_edge[b])/2
        axs30=axs3[0].plot(bins, hist/hist_norm,c=color[dd])








##### TOPPRED
#dd+=1
#get parameters name from the xml
mydoc = minidom.parse('bfm_sensitivity_toppred.xml')
items = mydoc.getElementsByTagName('parameter')
parameters = []
for i in range(len(items)):
    string = items[i].attributes['variable'].value
    if i<200:
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
    else :
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
#get targets names
target = []
items = mydoc.getElementsByTagName('target')
for i in range(len(items)):
    string = items[i].attributes['name'].value
    target.append(string)
#remove 05 03c 03h R3c
#### REMOVE ALSO O2O

key = 'O5'
keyOc= 'O3_c'
keyO2o= 'O2_o'
keyOh= 'O3h_h'
keyR3c = 'R3_c'
false_index = []
for i,name in enumerate(target):
    if name.find(key) !=-1 or name.find(keyOc) !=-1 or name.find(keyOh) !=-1 or name.find(keyR3c) !=-1 or name.find(keyO2o) !=-1 :    
                false_index.append(i)
#target.pop()

infile = open('bfm_sensitivity_toppred_extinction.pickle','rb')
new_dict = pickle.load(infile)
infile.close()
infile1 = open('bfm_sensitivity_toppred.pickle','rb')
new_dict1 = pickle.load(infile1)
infile1.close()

inputs = new_dict1.get('X')
outputs = new_dict.get('Y')
outputs = outputs.tolist()
#limit inputs to the folders generated
indexes = new_dict.get('I')
kk=0
for v in indexes:
    if v!=0:
        kk+=1
print('living: ',kk)
inputs = [inputs[int(i)] for i in indexes]
#outputs = [outputs[int(i)] for i in indexes]
outputs_living = np.copy(outputs).tolist()
target_living = np.copy(target).tolist()
for i in range(len(outputs)):  #remove O5
    for ind in sorted(false_index, reverse=True):
         outputs[i].pop(ind)
for ind in sorted(false_index, reverse=True):
    target.pop(ind)
count = np.zeros(len(inputs))
for o in range(len(outputs)):
    for i in range(int(len(outputs[o])/2)):
        if outputs[o][2*i]+outputs[o][2*i+1]==2:
    #if np.sum(outputs[o][1::2]) >1:
            count[o]=1
        else :
            pass#count[o]=0

#only living species
right_index_living = []
varnames = ["B1_c","P1_c","P2_c","P3_c","P4_c","Z5_c","Z3_c","Z4_c"]
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
print('non-stationary',lll)
print('toppred fraction',lll/len(outputs_living))

#histo only for living counts

count = count_living
outputs = outputs_living
target = target_living

x = np.zeros(len(inputs))
ind_n = []
for iv,var in enumerate(parameters):
    if var == 'N3n' or var=='N4n':
        ind_n.append(iv)

x_n=np.zeros(len(inputs))
for i in ind_n:
    for j in range(len(inputs)):
        x_n[j]+=inputs[j][i]
hist,bin_edge = np.histogram(x_n, weights=count)
hist_norm,bin_ed = np.histogram(x_n)
bins=np.zeros(len(hist))
for b in range(len(bin_edge)-1):
        bins[b]=(bin_edge[b+1]+bin_edge[b])/2
#axs31 = axs3[1].plot(bins, hist/hist_norm,c='blue',linestyle=dashList[dd],label='TOPPRED')
for i in range(len(inputs[0])):
    if parameters[i] == 'N1p':
        for j in range(len(inputs)):
            x[j]=inputs[j][i]
        hist,bin_edge = np.histogram(x, weights=count)
        hist_norm,bin_ed = np.histogram(x)
        bins=np.zeros(len(hist))
        for b in range(len(bin_edge)-1):
             bins[b]=(bin_edge[b+1]+bin_edge[b])/2
#       axs30=axs3[0].plot(bins, hist/hist_norm,c='blue',linestyle=dashList[dd])














for iax,ax in enumerate(fig3.get_axes()):
    ax.label_outer()
    ax.set_ylim(0,1. )
    my_xticks = ax.get_xticks()
    ax.set_yticks([0,0.5,1.])
    ax.axhline(y=0.5, color='r', linestyle='dotted')
    ax.set_xlabel('[$mmol m^{-3}$]')
    if iax == 0:
        ax.set_xticks([0.01,2.0])
        ax.annotate('N1p', xy=(0,1), xycoords='axes fraction', xytext=(1, -1), textcoords='offset points', horizontalalignment='left', verticalalignment='top')
    if iax == 1:
        ax.set_xticks([0.02,64.0])
        ax.annotate('N3n + N4n', xy=(0,1), xycoords='axes fraction', xytext=(1, -1), textcoords='offset points', horizontalalignment='left', verticalalignment='top')
fig3.legend(ncol=5, loc='upper center',fontsize=9)
fig3.savefig('living_PNhisto_color.png', format='png',dpi=350)
