#!/usr/bin/python

from fileio.fileread import *
import sys
import operator
import numpy as np
from iocontrol.options import get_options as getopt

options=['-i', '-o']
types=['str', 'str']
defaults=['file.in', 'default']

filename, outname=getopt(sys.argv, options, types, defaults)
if outname=='default':
    outname=filename.replace('.out', '_dev.out')


data, comment=file_to_array(filename)

x=data[:,0]
data=data[:,1:]


outdata=np.gradient(data, x, axis=0)
comment='#Derivatives w.r.t. first column'
outdata=np.c_[x,outdata]

np.savetxt(outname, outdata, fmt='%15.4f', delimiter='', newline='\n', header=comment, comments='')
print("Output is written in %s."%outname)

