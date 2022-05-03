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

temp=0.5
mass=1.0
rd.seed()
box=[0.0, 0.0, 0.0]
options=['-d', '-p', '-o', '-it', '-ot', '-v']
types=['str', 'str', 'str', 'str', 'str', 'float']
defaults=['drop.gro', 'brush.gro', 'conf.gro', 'topol_brush.top', 'topol.top', 2.0]
dropname, brushname, outname, intopol, outtopol, dropvel=getopt(sys.argv, options, types, defaults)


title_drop, topol_drop, box_drop, coord_drop, vel_drop=read_gro(dropname)
title_brush, topol_brush, box_brush, coord_brush, vel_brush=read_gro(brushname)
box_diff=[(x-y)/2 for x, y in zip(box_brush, box_drop)]
coord=[]
topol=[]
vel=[]
ndrop=len(coord_drop)
coord=[[coord_drop[i][0]+box_diff[0], coord_drop[i][1]+box_diff[1], coord_drop[i][2]+box_diff[2]] for i in range(ndrop)]
topol=topol_drop
vel=[[vel_drop[i][0], vel_drop[i][1], vel_drop[i][2]-dropvel] for i in range(ndrop)]
last_molnum=topol[-1][0]
last_atnum=topol[-1][3]
for i in range(len(coord_brush)):
    coord.append(coord_brush[i])
    topol.append([last_molnum+topol_brush[i][0], topol_brush[i][1], topol_brush[i][2], topol_brush[i][3]+last_atnum])
    vel.append(vel_brush[i])
title="Droplet over polymer brush"
write_gro(outname, title, topol, box_brush, coord, vel)
            
print(outname+" is generated.")

outfp=open(outtopol, "w")
infp=open(intopol, "r")
intopdata=infp.readlines()
idx=0
tobreak=False
while not tobreak and idx<len(intopdata):
    if intopdata[idx]=='[ System ]\n':
        outfp.write('[ System ]\n')
        idx+=1
        while idx<len(intopdata):
            if intopdata[idx]=='\n' or intopdata[idx][0]=='[':
                outfp.write('\n')
                tobreak=True
                break
            elif intopdata[idx][0]==';':
                outfp.write(intopdata[idx])
                outfp.write("%10s%10d\n"%("SOL", ndrop))
                idx+=1

            else:
                outfp.write(intopdata[idx])
                idx+=1
    else:
        idx+=1

outfp.write("[ SOL Atoms ]\n")
outfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
outfp.write("%10d%10s%10.1f\n"%(0, "W", 1.0))

idx=0
tobreak=False
while not tobreak and idx<len(intopdata):
    if intopdata[idx]=='[ POL Atoms ]\n':
        outfp.write(intopdata[idx])
        idx+=1
        while idx<len(intopdata):
            if intopdata[idx]=='\n' or intopdata[idx][0]=='[':
                outfp.write('\n')
                tobreak=True
                break
            else:
                outfp.write(intopdata[idx])
                idx+=1
    else:
        idx+=1

idx=0
tobreak=False
while not tobreak and idx<len(intopdata):
    if intopdata[idx]=='[ WAL Atoms ]\n':
        outfp.write(intopdata[idx])
        idx+=1
        while idx<len(intopdata):
            if intopdata[idx]=='\n' or intopdata[idx][0]=='[':
                outfp.write('\n')
                tobreak=True
                break
            else:
                outfp.write(intopdata[idx])
                idx+=1
    else:
        idx+=1


idx=0
tobreak=False
while not tobreak and idx<len(intopdata):
    if intopdata[idx]=='[ POL Bonds ]\n':
        outfp.write(intopdata[idx])
        idx+=1
        while idx<len(intopdata):
            if intopdata[idx]=='\n' or intopdata[idx][0]=='[':
                outfp.write('\n')
                tobreak=True
                break
            else:
                outfp.write(intopdata[idx])
                idx+=1
    else:
        idx+=1

idx=0
tobreak=False
while not tobreak and idx<len(intopdata):
    if intopdata[idx]=='[ Interbonds ]\n':
        outfp.write(intopdata[idx])
        idx+=1
        while idx<len(intopdata):
            if intopdata[idx]=='\n' or intopdata[idx][0]=='[':
                outfp.write('\n')
                tobreak=True
                break
            elif intopdata[idx][0]==';':
                outfp.write(intopdata[idx])
                idx+=1
            else:
                interbonds=intopdata[idx].split()
                outfp.write("%10d%10d%10s\n"%(int(interbonds[0])+ndrop, int(interbonds[1])+ndrop, interbonds[2]))
                idx+=1
    else:
        idx+=1

idx=0
tobreak=False
while not tobreak and idx<len(intopdata):
    if intopdata[idx]=='[ Nonbonded ]\n':
        outfp.write(intopdata[idx])
        idx+=1
        while idx<len(intopdata):
            if intopdata[idx]=='\n' or intopdata[idx][0]=='[':
                outfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("W", "W", 25, 0.75, -40, 1.0))
                outfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("W", "P", 25, 0.75, -40, 1.0))
                outfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("W", "S", 25, 0.75, -10, 1.0))
                outfp.write('\n')
                tobreak=True
                break
            else:
                outfp.write(intopdata[idx])
                idx+=1
    else:
        idx+=1

idx=0
tobreak=False
while not tobreak and idx<len(intopdata):
    if intopdata[idx]=='[ Bondlength ]\n':
        outfp.write(intopdata[idx])
        idx+=1
        while idx<len(intopdata):
            if intopdata[idx]=='\n' or intopdata[idx][0]=='[':
                outfp.write('\n')
                tobreak=True
                break
            else:
                outfp.write(intopdata[idx])
                idx+=1
    else:
        idx+=1



