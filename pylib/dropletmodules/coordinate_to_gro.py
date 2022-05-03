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
from dropletmodules.pbc import *




inname, outname=getopt(sys.argv, ['-i', '-o'], ['str', 'str'], ['coordinate', 'coord.gro'])
species, box, topol, coord, vel=read_coordinate(inname)


fpout=open(outname, "w")
nbeads=len(coord)
fpout.write("converted from coordinate\n")
fpout.write("%d\n"%nbeads)
for i in range(nbeads):
    if species[topol[i][3]-1][3]=='SOLID':
        name='SLD'
    elif species[topol[i][3]-1][3]=='SOLVENT':
        name='SVT'
    elif species[topol[i][3]-1][3]=='POLYMER':
        name='POL'
    else:
        name='NONE'

    fpout.write("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n"%(topol[i][1], name, 'T'+str(topol[i][3]), topol[i][2], get_particle_in_box(coord[i][0], box[0]), get_particle_in_box(coord[i][1], box[1]), get_particle_in_box(coord[i][2], box[2]), vel[i][0], vel[i][1], vel[i][2]))
fpout.write("%10.5f%10.5f%10.5f\n"%(box[0], box[1], box[2]))
fpout.close()
