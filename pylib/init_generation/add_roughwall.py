#!/usr/bin/python

#------------------------------------------------------------------------------#
#  A program to put solvent molecules into a box randomly producing an input of
# a MDPD program.
#------------------------------------------------------------------------------#

import sys
sys.path.append('../')
import numpy as np
from iocontrol.options import get_options as getopt
import random as rd
from coordinate import *
from pbc import *
from optparse import OptionParser



unitcell=0.5
options=['-b', '-w', '-i', '-o', '-x', '-y', '-z', '-gap', '-height', '-width', '-dir']
types=['float', 'float', 'str', 'str', 'float', 'float', 'float', 'float', 'float', 'float', 'float']
defaults=[10.0, 5.0, 'nowall.gro', 'wall.gro', 0, 0, 0, 0, 0, 0, 'xy']
nbox, wallgap, inname, outname, boxx, boxy, boxz, gap, height, width, direc=getopt(sys.argv, options, types, defaults)

title, topol, box, coord, vel=read_gro(inname)
newbox=[boxx, boxy, boxz+2*wallgap+height]
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
ptclidx=len(coord)+1
molidx=topol[-1][0]+1


numwall=0
nsp=int((gap+width)/unitcell)
gapunit=int(gap/unitcell)
basenzh=3
botbasez=wallgap-(basenzh-1)*unitcell
topbasez=wallgap+boxz+height+(basenzh-1)*unitcell

#You have to write the code to add rough pillared walls to the droplet only system
for i in range(nboxx):
    for j in range(nboxy):
        if nsp==0:
            nzheight=basenzh
        elif i%nsp<gapunit or j%nsp<gapunit:
            nzheight=basenzh
        else:
            nzheight=int(height/unitcell)+basenzh
        
        for k in range(nzheight):
            coord.append([i*unitcell, j*unitcell, botbasez+k*unitcell])
            vel.append([0,0,0])
            topol.append([molidx, "WAL", "W", ptclidx])
            ptclidx+=1
            molidx+=1
            numwall+=1
#Surface at top
for i in range(nboxx):
    for j in range(nboxy):
        if nsp==0:
            nzheight=basenzh
        elif i%nsp<gapunit or j%nsp<gapunit:
            nzheight==basenzh
        else:
            nzheight=int(height/unitcell)+basenzh

        for k in range(nzheight):
            coord.append([i*unitcell, j*unitcell, topbasez-k*unitcell])
            vel.append([0,0,0])
            topol.append([molidx, "WAL", "W", ptclidx])
            ptclidx+=1
            molidx+=1
            numwall+=1


title=title+" after rough wall added"
write_gro(outname, title, topol, newbox, coord, vel)
