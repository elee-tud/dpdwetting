#!/usr/bin/python

#------------------------------------------------------------------------------#
#  A program to put solvent molecules into a box randomly producing an input of
# a MDPD program.
#------------------------------------------------------------------------------#

import sys

sys.path.append("../")
import numpy as np
from iocontrol.options import get_options as getopt
import random as rd
from dropletmodules.coordinate import *
from dropletmodules.polymer_random_walk import *
from dropletmodules.pbc import *
from ring.ringconfiguration import *

temp=1.0
mass=1.0
rd.seed()
box=[0.0, 0.0, 0.0]
options=['-s', '-p', '-l', '-bl', '-b', '-o' , '-r', '-x', '-y', '-z', '-pp', '-ss', '-ps', '-ring']
types=['int', 'int', 'int', 'float', 'float', 'str', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'int']
defaults=[0, 0, 0, 1.0, 10.0, 'conf.gro', 0.2, 0., 0., 0., -40, -40, -40, 0]
nsol, npol, lpol, bl, boxl, outname, polgenboxratio, boxx, boxy, boxz, pp, ss, ps, isring=getopt(sys.argv, options, types, defaults)
if not boxx==0.:
    box=np.array([boxx, boxy, boxz])
else:
    box=np.array([boxl, boxl, boxl])
topol=[]
coord=[]
vel=[]
molidx=1
ptclidx=1
polgenbox=box*polgenboxratio
polgenmin=(box-polgenbox)/2

for i in range(npol):
    for j in range(lpol):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'POL', 'P', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    if isring==0:
        polcoord=rdpolymer(lpol, bl, polgenbox)
        for pc in polcoord:
            pc[0]+=polgenmin[0]
            pc[1]+=polgenmin[1]
            pc[2]+=polgenmin[2]
    else:
        polcoord=doublefolded(lpol, bl, box, polgenboxratio)
    coord+=polcoord

    molidx+=1
for i in range(len(coord)):
    coord[i]=get_particle_in_box(coord[i], box)
    
for i in range(nsol):
    topol.append([molidx, 'SOL', 'S', ptclidx])
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
if npol>0:
    topolfp.write("%10s%10d\n"%("POL", npol))
if nsol>0:
    topolfp.write("%10s%10d\n"%("SOL", nsol))

topolfp.write("\n")
topolfp.write(";Molecule Atom Information\n")
if nsol>0:
    topolfp.write("[ SOL Atoms ]\n")
    topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
    topolfp.write("%10d%10s%10.1f\n"%(0, "S", 1.0))
    topolfp.write("\n")
if npol>0:
    topolfp.write("[ POL Atoms ]\n")
    topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
    for i in range(lpol):
        topolfp.write("%10d%10s%10.1f\n"%(i, "P", 1.0))
    topolfp.write("\n")
    if lpol>1:
        topolfp.write(";Molecule Bond Information\n")
        topolfp.write("[ POL Bonds ]\n")
        topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
        for i in range(lpol-1):
            topolfp.write("%10d%10d%10s\n"%(i, i+1, "PP"))
        if not isring==0:
            topolfp.write("%10d%10d%10s\n"%(lpol-1, 0, "PP"))

        topolfp.write("\n")
topolfp.write(";Nonbonded Parameters\n")
topolfp.write("[ Nonbonded ]\n")
topolfp.write(";%9s%10s%10s%10s%10s%10s\n"%("Atype1","Atype2","B", "r_B", "A", "r_A"))
if nsol>0:
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "S", 25, 0.75, ss, 1.0))
    if npol>0:
        topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "P", 25, 0.75, ps, 1.0))
        topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("P", "P", 25, 0.75, pp, 1.0))
topolfp.write("\n")
if lpol>1:
    topolfp.write(";Bond length Parameters\n")
    topolfp.write("[ Bondlength ]\n")
    topolfp.write(";%9s%10s%10s\n"%("btype", "k", "l0"))
    topolfp.write("%10s%10.2f%10.3f\n"%("PP", 300.0, 0.65))
topolfp.close()


