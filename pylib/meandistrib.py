#!/home/elee/programs/Python-3.6.10/bin/python

from fileio.fileread import *
import sys
import operator
import numpy as np
from iocontrol.options import get_options as getopt

options=['-i']
types=['str']
defaults=['file.in']

[filename]=getopt(sys.argv, options, types, defaults)
print(filename)

data, comment=file_to_array(filename)

x=data[:,0]
data=data[:,1:]
ncol=np.size(data[0])
nrow=np.size(data[:,0])
avg=[0. for i in range(ncol)]
std=[0. for i in range(ncol)]
wsum=[0. for i in range(ncol)]
for i in range(ncol):
    for j in range(nrow):
        if not np.isnan(data[j,i]):
            wsum[i]=wsum[i]+data[j,i]
            avg[i]=avg[i]+x[j]*data[j,i]
            std[i]=std[i]+x[j]*x[j]*data[j,i]
    if not wsum[i]==0:
        avg[i]=avg[i]/wsum[i]
        std[i]=np.sqrt(std[i]/wsum[i]-avg[i]*avg[i])
    else:
        avg[i]=0
        std[i]=0
    print("Column %d-> mean: %15.5e  /  std: %15.5e"%(i+1, avg[i], std[i]))



