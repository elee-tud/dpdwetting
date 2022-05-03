#!/home/elee/programs/Python-3.6.10/bin/python

import sys
import numpy as np
from fileio.fileread import file_to_array as fta
from iocontrol.options import get_options




startname, outname, start, end=get_options(sys.argv, ['-i', '-o', '-s', '-e'], ['str', 'str', 'int', 'int'], ['file.in', 'averaged.xvg',  1, 1])
shape=[]
idx=0
data=[None for i in range(start, end+1)]
for i in range(start, end+1):
    filename=startname.replace(str(start), str(i), 1)
    print("processing "+filename+"...")
    
    data[idx], comment=fta(filename)
    shape.append(data[idx].shape)
    idx += 1
    
shape=np.array(shape)
alldata=[]
for d in data:
    alldata.append(d[:np.min(shape[:,0]),:])
alldata=np.array(alldata)
avg=np.mean(alldata, axis=0)
std=np.std(alldata, axis=0)
alldata=np.append(avg, std, 1)
alldata=np.delete(alldata, shape[0,1], 1)
comment='#Average of data in files '+startname.replace(str(start),'X')+"from X="+str(start)+" to "+str(end)+".\n"+comment





np.savetxt(outname, np.array(alldata), fmt='%15.4f', delimiter='', newline='\n', header=comment, comments='')
print("Data is written in %s."%outname)



