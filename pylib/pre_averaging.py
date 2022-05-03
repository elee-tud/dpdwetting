#!/home/elee/programs/Python-3.6.10/bin/python

from fileio.fileread import *
import sys
import operator
import numpy as np
from iocontrol.options import get_options as getopt

options=['-i', '-o', '-b', '-e', '-n']
types=['str', 'str', 'int', 'int', 'int']
defaults=['file.in', 'averaged.out',  0, -1, 5]

filename, outname, begstep, endstep, navg=getopt(sys.argv, options, types, defaults)

data, comment=file_to_array(filename)
if endstep==-1:
    endstep=len(data)

numcol=len(data[0])
numrow=len(data)
numoutrow=numrow-navg
out=[[0.0 for i in range(numcol)] for j in range(numoutrow)]
for i in range(numcol):
    for j in range(numoutrow):
        if j==0:
            for k in range(navg):
                out[j][i]+=data[j+k][i]
        else:
            out[j][i]=out[j-1][i]-data[j-1][i]+data[j+navg-1][i]

for i in range(numcol):
    for j in range(numoutrow):
        out[j][i]/=navg
comment='#Averaged data by N successive original data '

np.savetxt(outname, np.array(out), fmt='%15.4f', delimiter='', newline='\n', header=comment, comments='')
print("Output is written in %s."%outname)

