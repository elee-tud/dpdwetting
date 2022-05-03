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

temp=1.0
mass=1.0


rd.seed()

options=['-N', '-b', '-x', '-y', '-z', '-o'];
types=['int', 'float', 'float', 'float', 'float', 'str'];
defaults=[0, -1.0, 0.0, 0.0, 0.0, 'conf.gro']
box=[0,0,0]
nbeads, box_all, box[0], box[1], box[2], outname=getopt(sys.argv, options, types, defaults)

if box_all!=-1.0:
    box[0]=box_all
    box[1]=box_all
    box[2]=box_all

title='Pure solvent with '+str(nbeads)+' particles'
topol=[[i+1, 'SOL', 'W', i+1] for i in range(nbeads)]
coord=[[rd.random()*box[0], rd.random()*box[1], rd.random()*box[2]] for i in range(nbeads)]
vel=[np.random.normal(0.0, temp/mass, 3) for i in range(nbeads)]

write_gro(outname, title, topol, box, coord, vel)
print(outname+" is generated.")
