#!/home/elee/programs/Python-3.6.10/bin/python

from fileio.fileread import *
import sys
import operator
import numpy as np
from iocontrol.options import get_options as getopt

options=['-i', '-b', '-e', '-bz', '-dz', '-ez']
types=['str', 'int', 'int', 'float', 'float', 'float']
defaults=['file.in', 0, -1, 0.5, 1.0, 70.5]

filename, begstep, endstep, begz, dz, endz=getopt(sys.argv, options, types, defaults)

nz=int((endz-begz)/dz)+1
data, comment=file_to_array(filename)
z=[0.]
for i in range(nz):
    z.append(begz+dz*i)
data=np.append([z], data, axis=0)
newdata=data.transpose()
newdata=np.delete(newdata, 0, 0)
outname=filename.replace(".out", "_trans.out")
outheader="#Transposed"
np.savetxt(outname, newdata, fmt='%1.5e', newline='\n', header=outheader)

