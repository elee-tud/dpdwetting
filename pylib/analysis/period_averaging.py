#!/usr/bin/python

from fileio.fileread import *
import sys
import operator
import numpy as np
from iocontrol.options import get_options as getopt

options=['-i', '-o', '-b', '-e', '-p', '-n']
types=['str', 'str', 'int', 'int', 'float', 'int']
defaults=['file.in', 'period_averaged.out',  0, -1, 10, 1]

filename, outname, begstep, endstep, period, nbin=getopt(sys.argv, options, types, defaults)

data, comment=file_to_array(filename)
if endstep==-1:
    endstep=len(data)
data=np.array(data[begstep:endstep][:])
numx, numy=data.shape

gap=period/nbin
ndata=np.ndarray(nbin)

output=np.ndarray((nbin, numy))

for i in range(nbin):
    output[i,0]=gap*(i+0.5)

for idx, time in enumerate(data[:,0]):
    index=int((time%period)/gap)
    ndata[index]=ndata[index]+1
    for j in range(1,numy):
        output[index,j]=output[index, j]+data[idx, j]

for i in range(nbin):
    for j in range(1,numy):
        output[i][j]=output[i][j]/ndata[j]





comment='#Pre-Averaged data by with the period %f'%period

np.savetxt(outname, np.array(output), fmt='%15.4f', delimiter='', newline='\n', header=comment, comments='')
print("Output is written in %s."%outname)

