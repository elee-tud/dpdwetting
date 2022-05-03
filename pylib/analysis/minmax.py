#!/usr/bin/python

from fileio.fileread import *
import sys
import operator
import numpy as np
from iocontrol.options import get_options as getopt

options=['-i', '-b', '-e']
types=['str', 'int', 'int']
defaults=['file.in', 0, -1]

filename, begstep, endstep=getopt(sys.argv, options, types, defaults)

data, comment=file_to_array(filename)
if endstep==-1:
    endstep=len(data)
for i in range(1,len(data[0])):
    maxval=np.nanmax(data[begstep:endstep,i])
    maxndx=data[np.nanargmax(data[begstep:endstep,i]),0]
    minval=np.nanmin(data[begstep:endstep,i])
    minndx=data[np.nanargmin(data[begstep:endstep,i]),0]
    print("Column number %d: maximum( %8.5e , %8.5e ), minimum( %8.5e , %8.5e )"%(i+1, maxndx, maxval, minndx, minval))



