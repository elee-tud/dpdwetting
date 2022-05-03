#!/home/elee/programs/Python-3.6.10/bin/python

from fileio.fileread import *
import sys
import operator
import numpy as np
from iocontrol.options import get_options as getopt
from scipy.fft import fft, fftfreq

options=['-i', '-o', '-b', '-e']
types=['str', 'str', 'int', 'int']
defaults=['file.in', 'fft.out',  0, -1]

filename, outname, begstep, endstep=getopt(sys.argv, options, types, defaults)

data, comment=file_to_array(filename)
if endstep==-1:
    endstep=len(data)
data=np.array(data[begstep:endstep][:])
numx, numy=data.shape
print("number of x points=%d"%numx)

dt=data[1,0]-data[0,0]
print("time interval=%8.4e"%dt)

freq=fftfreq(numx, dt)
nfreq=numx//2

result=np.ndarray((nfreq-1, numy))
result[:,0]=freq[1:nfreq]
for j in range(1,numy):
    fftval=fft(data[:,j])
    result[:,j]=1/nfreq*np.abs(fftval[1:nfreq])




comment='#Furier transform'

np.savetxt(outname, np.array(result), fmt='%15.4f', delimiter='', newline='\n', header=comment, comments='')
print("Output is written in %s."%outname)

