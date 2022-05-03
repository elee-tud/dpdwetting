#!/usr/bin/python

from fileio.fileread import *
import sys
import operator
import numpy as np
from iocontrol.options import get_options as getopt

options=['-i', '-b', '-e', '-r']
types=['str', 'int', 'int', 'int']
defaults=['file.in', 0, -1, -1]

filename, begstep, endstep, rowlen=getopt(sys.argv, options, types, defaults)

data, comment=file_to_array(filename, rowlen)
if endstep==-1:
    endstep=len(data)
seldata=data[begstep:endstep, 1:]
avg=np.average(data[begstep:endstep], axis=0)
std=np.std(data[begstep:endstep], axis=0)

for i in range(1,len(avg)):
    print("Column number %d: average= %8.5e , std= %8.5e"%(i+1, avg[i], std[i]))
