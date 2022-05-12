import netCDF4 as nc
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('results')
args = parser.parse_args()

ds = nc.Dataset(args.results)
new = nc.Dataset(args.results)
substring = 'var'
for name in ds.variables.keys():
    if name.find(substring) != -1:
        break
    else : 
        new.variables.pop(name)

for name in new.variables.keys():
    print('<target name="' + name + '" expression="' + name + '" path="result.nc"/>')







