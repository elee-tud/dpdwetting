#!/home/elee/programs/Python-3.6.10/bin/python

#------------------------------------------------------------------------------#
#  A program to put solvent molecules into a box randomly producing an input of
# a MDPD program.
#------------------------------------------------------------------------------#

import sys
import numpy as np
from iocontrol.options import get_options as getopt
import random as rd
from dropletmodules.coordinate import *
from dropletmodules.polymer_random_walk import *
from dropletmodules.pbc import *

options=['-s', '-i', '-o', '-b', '-k', '-n']
types=['int', 'str', 'str', 'int', 'int', 'int']
defaults=[0, 'traj.gro', 'out.gro', 0, 0, 0]
step, inname, outname, beg, skip, numstep=getopt(sys.argv, options, types, defaults)

fpin=open(inname)
if numstep==0:
    for i in range(step-1):
        read_gro_att(fpin, skip=True)
    title, topol, box, coord, vel=read_gro_att(fpin, skip=False)
    write_gro(outname, title, topol, box, coord, vel)
else:
    final=beg+skip*(numstep-1)+1
    idx=0
    for i in range(final):
        title, topol, box, coord, vel=read_gro_att(fpin, skip=False)
        print(i)
        if i==beg+skip*idx:
            outname='conf'+str(idx+1)+'.gro'
            write_gro(outname, title, topol, box, coord, vel)
            idx+=1
            
            
            


