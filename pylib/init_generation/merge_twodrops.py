#!/home/elee/.bin/python

#------------------------------------------------------------------------------#
#  A program to merge two droplets into one simulation box
#------------------------------------------------------------------------------#

import sys
import numpy as np
from iocontrol.options import get_options as getopt
import random as rd
from dropletmodules.coordinate import *
from dropletmodules.pbc import *



wallgap=5.0
unitcell=0.5
options=['-i1', '-i2', '-o', '-v', '-s']
types=['str', 'str', 'str', 'float', 'float']
defaults=['drop1.gro', 'drop2.gro', 'conf.gro', 0.0, 1.0]
inname1, inname2, outname, impvel, spacing=getopt(sys.argv, options, types, defaults)

title1, topol1, box1, coord1, vel1=read_gro(inname1)
title2, topol2, box2, coord2, vel2=read_gro(inname2)

if not box1[1]==box2[1]:
    print ('Error:Box size along y-direction do not match.')
    exit()

if not box1[2]==box2[2]:
    print ('Error:Box size along z-direction do not match.')
    exit()

box=[box1[0]+box2[0], box1[1], box1[2]]

half_impvel=impvel/2
half_spacing=spacing/2

center1=[0., 0., 0.]
for c in coord1:
    center1=[center1[i]+c[i] for i in range(3)]
center1=[center1[i]/len(coord1) for i in range(3)]
cbox1=[box1[i]/2 for i in range(3)]
diff1=[center1[i]-cbox1[i] for i in range(3)]
print(center1, diff1)
ncoord1=[]
for idx, c in enumerate(coord1):
    ncoord1.append([(c[0]-diff1[0])%box1[0], (c[1]-diff1[1])%box1[1], (c[2]-diff1[2]+half_spacing)%box1[2]])


center2=[0., 0., 0.]
for c in coord2:
    center2=[center2[i]+c[i] for i in range(3)]
center2=[center2[i]/len(coord2) for i in range(3)]
cbox2=[box2[i]/2 for i in range(3)]
diff2=[center2[i]-cbox2[i] for i in range(3)]
print(center2, diff2)
ncoord2=[]
for idx, c in enumerate(coord2):
    ncoord2.append([(c[0]-diff2[0])%box2[0]+box1[0], (c[1]-diff2[1])%box2[1], (c[2]-diff2[2]-half_spacing)%box2[2]])

for idx in range(len(vel1)):
    vel1[idx][0]=vel1[idx][0]+half_impvel
for idx in range(len(vel2)):
    vel2[idx][0]=vel2[idx][0]-half_impvel

coord=ncoord1+ncoord2
vel=vel1+vel2


for i in range(len(topol2)):
    topol2[i][0]=topol2[i][0]+topol1[-1][0]
    topol2[i][3]=topol2[i][3]+topol1[-1][3]

topol=topol1+topol2
fcoord=[]
fvel=[]
ftopol=[]

for i in range(len(topol)):
    if topol[i][2].replace(" ", "")=="P":
        fcoord.append(coord[i])
        fvel.append(vel[i])
        ftopol.append(topol[i])

for i in range(len(topol)):
    if topol[i][2].replace(" ", "")=="S":
        fcoord.append(coord[i])
        fvel.append(vel[i])
        ftopol.append(topol[i])



title="Conformation of two dorplets"
write_gro(outname, title, ftopol, box, fcoord, fvel)
