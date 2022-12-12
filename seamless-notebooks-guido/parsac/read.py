import pickle
import glob
import numpy as np

#get paths
filenames = glob.glob('/g100_work/icei_Rosati/seamless-notebooks-guido/parsac/test_total/*.txt')
filenames.sort(key=lambda x: int(''.join(filter(str.isdigit, x)))) #sort by number

data = np.zeros(len(filenames))
nn = np.zeros(len(filenames))
for i,name in enumerate(filenames):
    with open(name, 'r') as f:
       helper  = np.array([float(line) for line in f])
       print(helper)
       try:
           helper = max(helper)
           data[i] = max(helper,0)
       except:
           data[i]=-1
    nn[i] = int((name.rsplit('/')[-1]).rsplit('.')[-2])

#pickling
pkname = 'lyapunov_total.pickle'
outfile = open(pkname,'wb')
out_dict = {'L': data, 'N':nn}
pickle.dump(out_dict,outfile)
outfile.close()

chaos = 0
err = 0
read = 0
for ll in data:
    if ll == 0:
        read+=1
    if ll>0.001:
        chaos+=1
    if ll == -1:
        err+=1
print('chaos',chaos)
print('err',err)
print('done',read)
