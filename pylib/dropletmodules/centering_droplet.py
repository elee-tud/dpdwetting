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



options=['-i', '-o', '-b']
types=['str', 'str', 'float']
defaults=['before.gro', 'centered.gro', -1.0]
inname, outname, nbox=getopt(sys.argv, options, types, defaults)
title, topol, box, coord, vel=read_gro(inname)
if nbox==-1.0:
    newbox=box
else:
    newbox=np.array([nbox, nbox, nbox])

refrange=[[box[0]/2-1.0, box[0]/2+1.0],[box[1]/2-1.0, box[1]/2+1.0],[box[2]/2-1.0, box[2]/2+1.0]]
ref=np.array(coord[1])
for i in range(len(coord)):
    if coord[i][0]>=refrange[0][0] and coord[i][0]<refrange[0][1] and coord[i][1]>=refrange[1][0] and coord[i][1]<refrange[1][1] and coord[i][2]>=refrange[2][0] and coord[i][2]<refrange[2][1]:
        ref=np.array(coord[i])
        break

com=np.array([0.0,0.0,0.0])
for i in range(len(coord)):
    rij=get_shortest_image_vector(np.array(coord[i])-ref, box)
    com=com+rij


    

com=com/len(coord)+ref

for i in range(len(coord)):
    coord[i]=get_particle_in_box(np.array(coord[i])-com+np.array(newbox)/2, np.array(newbox))

title=title+" after centered"
write_gro(outname, title, topol, newbox, coord, vel)
