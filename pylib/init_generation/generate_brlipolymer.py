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
from ring.ringconfiguration import *

temp=1.0
mass=1.0
rd.seed()
box=[0.0, 0.0, 0.0]
options=['-na', '-la', '-nb', '-lb', '-nbr', '-lbr', '-nbrp', '-lbrp', '-bl', '-o' , '-x', '-y', '-z']
types=['int', 'int', 'int', 'int', 'int', 'int', 'int', 'int','float', 'str', 'float', 'float', 'float']
defaults=[0, 0, 0, 0, 0, 0, 0, 0, 0.6, 'conf.gro', 10.0, 10.0, 10.0]
npa, lpa, npb, lpb, npc, lpc, nbrc, lbrc, blength, outname, boxx, boxy, boxz=getopt(sys.argv, options, types, defaults)
box=np.array([boxx, boxy, boxz])
topol=[]
coord=[]
vel=[]
molidx=1
ptclidx=1

dbr=int(lpc/nbrc)

totna=npa*lpa
totnb=npb*lpb
totnc=npc*(lpc+nbrc*lbrc)
ntotbeads=totna+totnb+totnc
halfnc=int(npc/2)

abox=float(totna/ntotbeads)*box[0]
brbox=totnc/2/ntotbeads*box[0]
bbox=totnb/ntotbeads*box[0]


polcoord=[]
genbox=[abox, box[1], box[2]]
for i in range(npa):
    for j in range(lpa):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'APOL', 'A', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    polcoord+=rdpolymer(lpa, blength, genbox)
    molidx+=1
for i in range(len(polcoord)):
    coord.append(get_particle_in_box([polcoord[i][0]-abox/2, polcoord[i][1], polcoord[i][2]], box))

polcoord=[]
genbox=[bbox, box[1], box[2]]
for i in range(npb):
    for j in range(lpb):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'BPOL', 'B', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    polcoord+=rdpolymer(lpa, blength, genbox)
    molidx+=1
for i in range(len(polcoord)):
    coord.append(get_particle_in_box([polcoord[i][0]+abox/2+brbox, polcoord[i][1], polcoord[i][2]], box))

polcoord=[]
genbox=[brbox, box[1], box[2]]
for i in range(halfnc):
    for j in range(lpc):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'BRAN', 'A', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    for j in range(nbrc*lbrc):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'BRAN', 'B', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1

    polcoord+=branchedpolymer(lpc, dbr, lbrc, blength, genbox)
    molidx+=1

for i in range(len(polcoord)):
    coord.append(get_particle_in_box([polcoord[i][0]+abox/2, polcoord[i][1], polcoord[i][2]], box))


polcoord=[]
genbox=[brbox, box[1], box[2]]
for i in range(halfnc):
    for j in range(lpc):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'BRAN', 'A', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    for j in range(nbrc*lbrc):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'BRAN', 'B', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    polcoord+=branchedpolymer(lpc, dbr, lbrc, blength, genbox)

    molidx+=1

for i in range(len(polcoord)):
    coord.append(get_particle_in_box([polcoord[i][0]+abox/2+brbox+bbox, polcoord[i][1], polcoord[i][2]], box))

print(abox/2, abox/2+brbox, abox/2+brbox+bbox)

title="Polymer with solvent"
write_gro(outname, title, topol, box, coord, vel)
print(outname+" is generated.")


topolfp=open("topol.top", "w")
topolfp.write(";System information\n")
topolfp.write("[ System ]\n")
topolfp.write(";%9s%10s\n"%("Mol_name","Num_of_mols"))
topolfp.write("%10s%10d\n"%("APOL", npa))
topolfp.write("%10s%10d\n"%("BPOL", npb))
topolfp.write("%10s%10d\n"%("BRAN", npc))
topolfp.write("\n")
topolfp.write(";Molecule Atom Information\n")
topolfp.write("\n")

topolfp.write("[ APOL Atoms ]\n")
topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
for i in range(lpa):
    topolfp.write("%10d%10s%10.1f\n"%(i, "A", 1.0))
topolfp.write("\n")

topolfp.write(";Molecule Bond Information\n")
topolfp.write("[ APOL Bonds ]\n")
topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
for i in range(lpa-1):
    topolfp.write("%10d%10d%10s\n"%(i, i+1, "PP"))
topolfp.write("\n")

topolfp.write("[ BPOL Atoms ]\n")
topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
for i in range(lpb):
    topolfp.write("%10d%10s%10.1f\n"%(i, "B", 1.0))
topolfp.write("\n")

topolfp.write(";Molecule Bond Information\n")
topolfp.write("[ BPOL Bonds ]\n")
topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
for i in range(lpb-1):
    topolfp.write("%10d%10d%10s\n"%(i, i+1, "PP"))
topolfp.write("\n")

topolfp.write("[ BRAN Atoms ]\n")
topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
for i in range(lpc):
    topolfp.write("%10d%10s%10.1f\n"%(i, "A", 1.0))
for i in range(nbrc*lbrc):
    topolfp.write("%10d%10s%10.1f\n"%(i+lpc, "B", 1.0))


topolfp.write(";Molecule Bond Information\n")
topolfp.write("[ BRAN Bonds ]\n")
topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
for i in range(lpc-1):
    topolfp.write("%10d%10d%10s\n"%(i, i+1, "PP"))
for i in range(nbrc):
    topolfp.write("%10d%10d%10s\n"%((i+1)*dbr-1, lpc+i*lbrc, "PP"))
    for j in range(lbrc-1):
        topolfp.write("%10d%10d%10s\n"%(lpc+i*lbrc+j, lpc+i*lbrc+j+1, "PP"))

topolfp.write("\n")
topolfp.write(";Nonbonded Parameters\n")
topolfp.write("[ Nonbonded ]\n")
topolfp.write(";%9s%10s%10s%10s\n"%("Atype1","Atype2","B", "r_B"))
topolfp.write("%10s%10s%10.2f%10.2f\n"%("A", "A", 25, 1.0))
topolfp.write("%10s%10s%10.2f%10.2f\n"%("B", "B", 25, 1.0))
topolfp.write("%10s%10s%10.2f%10.2f\n"%("A", "B", 40, 1.0))
topolfp.write("\n")

topolfp.write(";Bond length Parameters\n")
topolfp.write("[ Bondlength ]\n")
topolfp.write(";%9s%10s%10s\n"%("btype", "k", "l0"))
topolfp.write("%10s%10.2f%10.3f\n"%("PP", 2.0, 0.0))
topolfp.close()

