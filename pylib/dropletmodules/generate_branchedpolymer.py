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
from ring.ringconfiguration import *

temp=1.0
mass=1.0
rd.seed()
box=[0.0, 0.0, 0.0]
options=['-p', '-l', '-db', '-lb', '-bl', '-o' , '-x', '-y', '-z']
types=['int', 'int', 'int', 'int', 'float', 'str', 'float', 'float', 'float']
defaults=[0, 0, 0, 0, 0.6, 'conf.gro', 10.0, 10.0, 10.0]
npol, lbb, dbr, lbr, blength, outname, boxx, boxy, boxz=getopt(sys.argv, options, types, defaults)
box=np.array([boxx, boxy, boxz])
topol=[]
coord=[]
vel=[]
molidx=1
ptclidx=1

nbr=int(lbb/dbr)
nbeads=lbb+nbr*lbr
for i in range(npol):
    for j in range(nbeads):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'POL', 'P', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    coord+=branchedpolymer(lbb, dbr, lbr, blength, box)


    molidx+=1
for i in range(len(coord)):
    coord[i]=get_particle_in_box(coord[i], box)
    

title="Polymer with solvent"
write_gro(outname, title, topol, box, coord, vel)
print(outname+" is generated.")

topolfp=open("topol.top", "w")
topolfp.write(";System information\n")
topolfp.write("[ System ]\n")
topolfp.write(";%9s%10s\n"%("Mol_name","Num_of_mols"))
topolfp.write("%10s%10d\n"%("POL", npol))
topolfp.write("\n")
topolfp.write(";Molecule Atom Information\n")
topolfp.write("\n")
topolfp.write("[ POL Atoms ]\n")
topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
for i in range(nbeads):
    topolfp.write("%10d%10s%10.1f\n"%(i, "P", 1.0))
topolfp.write("\n")
topolfp.write(";Molecule Bond Information\n")
topolfp.write("[ POL Bonds ]\n")
topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
for i in range(lbb-1):
    topolfp.write("%10d%10d%10s\n"%(i, i+1, "PP"))
for i in range(nbr):
    topolfp.write("%10d%10d%10s\n"%((i+1)*dbr-1, lbb+i*lbr, "PP"))
    for j in range(lbr-1):
        topolfp.write("%10d%10d%10s\n"%(lbb+i*lbr+j, lbb+i*lbr+j+1, "PP"))

topolfp.write("\n")
topolfp.write(";Nonbonded Parameters\n")
topolfp.write("[ Nonbonded ]\n")
topolfp.write(";%9s%10s%10s%10s%10s%10s\n"%("Atype1","Atype2","B", "r_B", "A", "r_A"))
topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("P", "P", 25, 1.0, 0.0, 0.0))
topolfp.write("\n")
topolfp.write(";Bond length Parameters\n")
topolfp.write("[ Bondlength ]\n")
topolfp.write(";%9s%10s%10s\n"%("btype", "k", "l0"))
topolfp.write("%10s%10.2f%10.3f\n"%("PP", 2.0, 0.0))
topolfp.close()


