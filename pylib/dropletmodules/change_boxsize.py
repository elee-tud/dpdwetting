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



wallgap=5.0
unitcell=1.0
nbox, inname, outname=getopt(sys.argv, ['-b', '-i', '-o'], ['float', 'str', 'str'], [10.0, 'old_coord', 'coordinate'])
species, box, topol, coord, vel=read_coordinate(inname)

newbox=[nbox, nbox, nbox]
diffx=(newbox[0]-box[0])/2
diffy=(newbox[1]-box[1])/2
diffz=(newbox[2]-box[2])/2
nboxx=int(newbox[0]/unitcell)
nboxy=int(newbox[1]/unitcell)
nbeads_wall=2*nboxx*nboxy
for i in range(len(coord)):
    coord[i][0]=get_particle_in_box(coord[i][0], box[0])+diffx
    coord[i][1]=get_particle_in_box(coord[i][1], box[1])+diffy
    coord[i][2]=get_particle_in_box(coord[i][2], box[2])+diffz


write_coordinate(outname, species, newbox, topol, coord, vel)
