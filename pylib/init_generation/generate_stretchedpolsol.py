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
from dropletmodules.polymer_random_walk import *
from dropletmodules.pbc import *
from pbc import *

temp=0.5
mass=1.0
rd.seed()
box=[0.0, 0.0, 0.0]
options=['-s', '-p', '-l', '-bl', '-b', '-o' ]
types=['int', 'int', 'int', 'float', 'float', 'str']
defaults=[0, 0, 0, 0.65, 10.0, 'conf.gro']
nsol, npol, lpol, bl, boxx, outname=getopt(sys.argv, options, types, defaults)
box=np.array([boxx, boxx, boxx])
topol=[]
coord=[]
vel=[]
molidx=1
ptclidx=1

firstpol=np.ndarray(shape=(npol, 3), dtype=float)

for i in range(npol):
    for j in range(lpol):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'POL', 'P', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    polcoord=np.ndarray(shape=(lpol, 3), dtype=float)
    while True:
        first=np.array([rd.random()*box[0], rd.random()*box[1], rd.random()*box[2]])
        overlap=False
        for k in range(j-1):
            vec=get_shortest_image_vector(first-firstpol[k], box)
            if (vec[0]*vec[0]+vec[1]*vec[1])**0.5 < 0.65:
                overlap=True
                break
        if not overlap:
            polcoord[0]=first
            firstpol[i]=first
            break
    for j in range(1,lpol):
        polcoord[j]=polcoord[j-1]+np.array([bl, 0, 0])
    coord+=polcoord.tolist()
    print("Polymer %d is generated"% (i+1))

    molidx+=1
for i in range(len(coord)):
    coord[i]=get_particle_in_box(coord[i], box)
    
for i in range(nsol):
    topol.append([molidx, 'SOL', 'W', ptclidx])
    coord.append([rd.random()*box[0], rd.random()*box[1],rd.random()*box[2]])
    vel.append(np.random.normal(0.0, temp/mass, 3))
    ptclidx+=1
    molidx+=1

title="Polymer with solvent"
write_gro(outname, title, topol, box, coord, vel)
print(outname+" is generated.")

topolfp=open("topol.top", "w")
topolfp.write(";System information\n")
topolfp.write("[ System ]\n")
topolfp.write(";%9s%10s\n"%("Mol_name","Num_of_mols"))
topolfp.write("%10s%10d\n"%("POL", npol))
if nsol>0:
    topolfp.write("%10s%10d\n"%("SOL", nsol))

topolfp.write("\n")
topolfp.write(";Molecule Atom Information\n")
if nsol>0:
    topolfp.write("[ SOL Atoms ]\n")
    topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
    topolfp.write("%10d%10s%10.1f\n"%(0, "W", 1.0))
    topolfp.write("\n")
topolfp.write("[ POL Atoms ]\n")
topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
for i in range(lpol):
    topolfp.write("%10d%10s%10.1f\n"%(i, "P", 1.0))
topolfp.write("\n")
topolfp.write(";Molecule Bond Information\n")
topolfp.write("[ POL Bonds ]\n")
topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
for i in range(lpol-1):
    topolfp.write("%10d%10d%10s\n"%(i, i+1, "PP"))
topolfp.write("\n")
topolfp.write(";Nonbonded Parameters\n")
topolfp.write("[ Nonbonded ]\n")
topolfp.write(";%9s%10s%10s%10s%10s%10s\n"%("Atype1","Atype2","B", "r_B", "A", "r_A"))
if nsol>0:
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("W", "W", 25, 0.75, -40, 1.0))
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("W", "P", 25, 0.75, -40, 1.0))
topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("P", "P", 25, 0.75, -40, 1.0))
topolfp.write("\n")
if lpol>1:
    topolfp.write(";Bond length Parameters\n")
    topolfp.write("[ Bondlength ]\n")
    topolfp.write(";%9s%10s%10s\n"%("btype", "k", "l0"))
    topolfp.write("%10s%10.2f%10.3f\n"%("PP", 300.0, 0.65))
topolfp.close()


