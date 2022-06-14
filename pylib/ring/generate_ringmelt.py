#!/p/home/jusers/lee8/juwels/bin/python

#------------------------------------------------------------------------------#
#  A program to put solvent molecules into a box randomly producing an input of
# a MDPD program.
#------------------------------------------------------------------------------#

import sys
import numpy as np
from iocontrol.options import get_options as getopt
import random as rd
from dropletmodules.coordinate import *
from ring.ringconfiguration import *
from dropletmodules.pbc import *

temp=0.5
mass=1.0
rd.seed()
box=[0.0, 0.0, 0.0]
options=['-p', '-l', '-bl', '-b', '-o', '-x', '-y', '-z']
types=['int', 'int', 'float', 'float', 'str', 'float', 'float', 'float']
defaults=[0, 0, 1.0, 10.0, 'conf.gro', 0., 0., 0.]
npol, lpol, bl, boxl, outname,  boxx, boxy, boxz=getopt(sys.argv, options, types, defaults)
if not boxx==0.:
    box=np.array([boxx, boxy, boxz])
else:
    box=np.array([boxl, boxl, boxl])
topol=[]
coord=[]
vel=[]
molidx=1
ptclidx=1
maxinteg=int(box[0]*1000)+1
occupied=[[False, False, False] for i in range(maxinteg)]


for i in range(npol):
    for j in range(lpol):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'POL', 'P', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    polcoord=doublefolded(lpol, bl, box)

    coord+=polcoord

    molidx+=1
for i in range(len(coord)):
    coord[i]=get_particle_in_box(coord[i], box)
    

title="Ring polymer melt"
write_gro(outname, title, topol, box, coord, vel)
print(outname+" is generated.")

topolfp=open("topol.top", "w")
topolfp.write(";System information\n")
topolfp.write("[ System ]\n")
topolfp.write(";%9s%10s\n"%("Mol_name","Num_of_mols"))
topolfp.write("%10s%10d\n"%("POL", npol))

topolfp.write("\n")
topolfp.write(";Molecule Atom Information\n")
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
topolfp.write("%10d%10d%10s\n"%(lpol-1, 0, "PP"))
topolfp.write("\n")
topolfp.write(";Nonbonded Parameters\n")
topolfp.write("[ Nonbonded ]\n")
topolfp.write(";%9s%10s%10s%10s\n"%("Atype1","Atype2","B", "r_B"))
topolfp.write("%10s%10s%10.2f%10.2f\n"%("P", "P", 25, 1.0))
topolfp.write("\n")
if lpol>1:
    topolfp.write(";Bond length Parameters\n")
    topolfp.write("[ Bondlength ]\n")
    topolfp.write(";%9s%10s%10s\n"%("btype", "k", "l0"))
    topolfp.write("%10s%10.2f%10.3f\n"%("PP", 2.0, 0.0))
topolfp.close()


