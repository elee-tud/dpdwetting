#!/usr/bin/python

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
unitcell=0.5
options=['-b', '-w', '-i', '-o', '-x', '-y', '-z', '-gap', '-height', '-width']
types=['float', 'float', 'str', 'str', 'float', 'float', 'float', 'float', 'float', 'float']
defaults=[10.0, 5.0, 'nowall.gro', 'wall.gro', 0, 0, 0, 0, 0, 0]
nbox, wallgap, inname, outname, boxx, boxy, boxz, gap, height, width=getopt(sys.argv, options, types, defaults)

title, topol, box, coord, vel=read_gro(inname)
newbox=[nbox, nbox, nbox+2*wallgap]
if boxx!=0:
    newbox[0]=boxx
if boxy!=0:
    newbox[1]=boxy
if boxz!=0:
    newbox[2]=boxz+2*wallgap
diffx=(newbox[0]-box[0])/2
diffy=(newbox[1]-box[1])/2
diffz=(newbox[2]-box[2])/2
nboxx=int(newbox[0]/unitcell)
nboxy=int(newbox[1]/unitcell)
for i in range(len(coord)):
    coord[i][0]=get_particle_vector_in_box(coord[i][0], box[0])+diffx
    coord[i][1]=get_particle_vector_in_box(coord[i][1], box[1])+diffy
    coord[i][2]=get_particle_vector_in_box(coord[i][2], box[2])+diffz
last_molnum=topol[-1][0]
last_atnum=topol[-1][0]
index=1

nxsp=int((gap+width)/unitcell)
nysp=int((gap+width)/unitcell)
gapunit=int(gap/unitcell)

for i in range(nboxx):
    for j in range(nboxy):
        if (int(i%nxsp)<gapunit or (j%nysp)<gapunit):
            nzht=2
        else:
            nzht=int(height/unitcell)+2
        for k in range(nzht):
            coord.append([i*unitcell, j*unitcell, wallgap-unitcell+k*unitcell])
            vel.append([0.0, 0.0, 0.0])
            topol.append([last_molnum+index, "WAL", "S", last_atnum+index])
            index+=1


for i in range(nboxx):
    for j in range(nboxy):
        coord.append([i*unitcell, j*unitcell, newbox[2]-wallgap])
        vel.append([0.0, 0.0, 0.0])
        topol.append([last_molnum+index, "WAL", "S", last_atnum+index])
        index+=1
for i in range(nboxx):
    for j in range(nboxy):
        coord.append([i*unitcell, j*unitcell, newbox[2]-wallgap+unitcell])
        vel.append([0.0, 0.0, 0.0])
        topol.append([last_molnum+index, "WAL", "S", last_atnum+index])
        index+=1

nbeads_wall=index

title=title+" after rough wall added"
write_gro(outname, title, topol, newbox, coord, vel)
