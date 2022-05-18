import pickle
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import argparse
from xml.dom import minidom
import glob

parser = argparse.ArgumentParser()
parser.add_argument('pickle')
parser.add_argument('xml')
args = parser.parse_args()

#get parameters name from the xml
mydoc = minidom.parse(args.xml)
items = mydoc.getElementsByTagName('parameter')
parameters = [] 
for i in range(len(items)):
    string = items[i].attributes['variable'].value
    if i<200:
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
    else :
        parameters.append(string.rsplit('/')[1]+string.rsplit('/')[-1])
#filenames = glob.glob('/g100_scratch/userexternal/gocchipi/BFM_TOTAL/*/*.nc')
#filenames.sort(key=lambda x: int(''.join(filter(str.isdigit, x)))) #sort by number

#indexes = np.zeros(len(filenames))
#for iv,var in enumerate(filenames):
#     indexes[iv] = int(var.rsplit('/')[-2])

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
    if name.find(key) !=-1 or name.find(keyOc) !=-1 or name.find(keyOh) !=-1 or name.find(keyR3c) !=-1 or name.find(keyO2o):
        false_index.append(i)

infile = open(args.pickle,'rb')
new_dict = pickle.load(infile)
infile.close()
#infile1 = open('bfm_sensitivity_total.pickle','rb')
#new_dict1 = pickle.load(infile1)
#infile1.close()
fig, axs = plt.subplots(10, 10)
fig1, ax1 = plt.subplots(10,10)
fig2, ax2 = plt.subplots(6,10)
plt.rcParams.update({'font.size': 4})
col = 0
col1 = 0
col2 = 0
inputs = new_dict.get('X')
outputs = new_dict.get('Y')
outputs = outputs.tolist()
#limit inputs to the folders generated
#indexes = new_dict.get('I')
#inputs = [inputs[int(i)] for i in indexes]


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
x = np.zeros(len(inputs))
for i in range(len(inputs[0])):
    for j in range(len(inputs)):
        x[j]=inputs[j][i]
    hist,bin_edge = np.histogram(x, weights=count)
    hist_norm,bin_ed = np.histogram(x)
    bins=np.zeros(len(hist))
    for b in range(len(bin_edge)-1):
        bins[b]=(bin_edge[b+1]+bin_edge[b])/2
    if i < 10 :
        ax = axs[col, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 20 :
        ax = axs[col+1, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 30 :
        ax = axs[col+2, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 40 :
        ax =axs[col+3, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 50 :
        ax = axs[col+4, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 60 :
        ax = axs[col+5, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 70 :
        ax = axs[col+6, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 80 :
        ax = axs[col+7, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 90 :
        ax = axs[col+8, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 100 :
        ax = axs[col+9, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 110 :
        ax = ax1[col1, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 120 :
        ax = ax1[col1+1, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 130 :
        ax = ax1[col1+2, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 140 :
        ax = ax1[col1+3, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 150 :
        ax = ax1[col1+4, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 160 :
        ax = ax1[col1+5, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 170 :
        ax = ax1[col1+6, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 180 :
        ax = ax1[col1+7, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 190 :
        ax = ax1[col1+8, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 200 :
        ax = ax1[col1+9, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 210 :
        ax = ax2[col2, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 220 :
        ax = ax2[col2+1, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 230 :
        ax = ax2[col2+2, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 240 :
        ax = ax2[col2+3, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 250 :
        ax = ax2[col2+4, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 260 :
        ax = ax2[col2+5, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 270 :
        ax = ax2[col2+6, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 280 :
        ax = ax2[col2+7, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 290 :
        ax = ax2[col2+8, i%10].plot(bins, hist/hist_norm, c='black')
    elif i < 300 :
        ax = ax2[col2+9, i%10].plot(bins, hist/hist_norm, c='black')

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

ax2[5,4].plot(bins, hist/hist_norm,c='black')

for iax,ax in enumerate(fig.get_axes()):
    ax.label_outer()
    ax.axhline(y=0.05, color='r', linestyle="dotted")
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    if iax < len(parameters):
        ax.annotate(parameters[iax], xy=(0,1), xycoords='axes fraction', xytext=(1, -1), textcoords='offset points', horizontalalignment='left', verticalalignment='top')
for iax,ax in enumerate(fig1.get_axes()):
    iax1 = iax + len(fig.get_axes())
    ax.label_outer()
    ax.axhline(y=0.5, color='r', linestyle='dotted')
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.annotate(parameters[iax1], xy=(0,1), xycoords='axes fraction', xytext=(1, -1), textcoords='offset points', horizontalalignment='left', verticalalignment='top')
for iax,ax in enumerate(fig2.get_axes()):
    iax1 = iax + len(fig.get_axes()) + len(fig1.get_axes())
    ax.label_outer()
    ax.set_ylim(0,1 )
    ax.set_xticks([])
    ax.set_yticks([])
    if iax1 < len(parameters):
        ax.annotate(parameters[iax1], xy=(0,1), xycoords='axes fraction', xytext=(1, -1), textcoords='offset points', horizontalalignment='left', verticalalignment='top')
        ax.axhline(y=0.5, color='r', linestyle='dotted')
    elif iax1 == len(parameters):
        ax.annotate('N3n + N4n', xy=(0,1), xycoords='axes fraction', xytext=(1, -1), textcoords='offset points', horizontalalignment='left', verticalalignment='top')
        ax.axhline(y=0.5, color='r', linestyle='dotted')


fig3, axs3 = plt.subplots()
axs3.plot(bins, hist/hist_norm,c='black')
axs3.axhline(y=0.5, color='r', linestyle="dotted")
string_name = parameters[ind_n[0]]+' + '+parameters[ind_n[1]]
axs3.set_title('N3n + N4n')

fig.savefig('FINAL_total001-100_ymax0dot1.png', format='png',dpi=150)
fig1.savefig('FINAL_total101-200_ymax0dot1.png', format='png',dpi=150)
fig2.savefig('FINAL_total201-300_ymax0dot1.png', format='png',dpi=150)


fig3.savefig('N3nN4n.png', format='png', dpi=150)
