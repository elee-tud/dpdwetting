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
options=['-l', '-bl', '-d', '-o' , '-x', '-y', '-z', '-w', '-u']
types=['int', 'float', 'float', 'str', 'float', 'float', 'float', 'float', 'float']
defaults=[0, 0.65, 0.0, 'conf.gro',  0., 0., 0., 5.0, 0.5]
lpol, bl, gden,  outname, boxx, boxy, boxz, wallgap, unitcell=getopt(sys.argv, options, types, defaults)


box=np.array([boxx, boxy, boxz])
topol=[]
coord=[]
vel=[]
molidx=1
ptclidx=1
nboxx=int(box[0]/unitcell)
nboxy=int(box[1]/unitcell)
wallcoord=[]
for i in range(nboxx):
    for j in range(nboxy):
        wallcoord.append([i*unitcell, j*unitcell, wallgap])
for i in range(nboxx):
    for j in range(nboxy):
        wallcoord.append([i*unitcell, j*unitcell, boxz-wallgap])

nwall=len(wallcoord)
npol=int(gden*box[0]*box[1])

occupied=[]
boundary=wallgap+bl
for i in range(npol):
    while True:
        thp=[int(rd.random()*nboxx), int(rd.random()*nboxy)]
        if thp not in occupied:
            polcoord=rdpolymer(lpol, bl, box, first=[thp[0]*unitcell+rd.uniform(-0.1,0.1), thp[1]*unitcell+rd.uniform(-0.1,0.1), boundary])
            occupied.append(thp)
            for j in range(lpol):
                if polcoord[j][2]<boundary:
                    polcoord[j][2]=2*boundary-polcoord[j][2]
                polcoord[j]=get_particle_in_box(polcoord[j], box)
            break
    for j in range(lpol):
        coord.append(polcoord[j])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'POL', 'P', ptclidx])
        ptclidx+=1
    molidx+=1
    

for i in range(len(wallcoord)):
    coord.append(wallcoord[i])
    vel.append([0.,0.,0.])
    topol.append([molidx, "WAL", "S", ptclidx])
    ptclidx+=1
    molidx+=1

            

    

title="polymer grafting wall"
write_gro(outname, title, topol, box, coord, vel)
print(outname+" is generated.")

topolfp=open("topol.top", "w")
topolfp.write(";System information\n")
topolfp.write("[ System ]\n")
topolfp.write(";%9s%10s\n"%("Mol_name","Num_of_mols"))
topolfp.write("%10s%10d\n"%("POL", npol))
topolfp.write("%10s%10d%10s\n"%("WAL", nwall, "frozen"))

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
topolfp.write("\n")

topolfp.write("[ WAL Atoms ]\n")
topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
topolfp.write("%10d%10s%10.1f\n"%(0, "S", 1.0))
topolfp.write("\n")

topolfp.write(";Nonbonded Parameters\n")
topolfp.write("[ Nonbonded ]\n")
topolfp.write(";%9s%10s%10s%10s%10s%10s\n"%("Atype1","Atype2","B", "r_B", "A", "r_A"))
topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "S", 25, 0.75, -40, 1.0))
topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "P", 25, 0.75, -10, 1.0))
topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("P", "P", 25, 0.75, -40, 1.0))
topolfp.write("\n")
topolfp.write("[ Interbonds ]\n")
topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
for i in range(npol):
    topolfp.write("%10d%10d%10s\n"%(i*lpol, npol*lpol+occupied[i][0]*nboxy+occupied[i][1], "PP"))
topolfp.write("\n")

topolfp.write(";Bond length Parameters\n")
topolfp.write("[ Bondlength ]\n")
topolfp.write(";%9s%10s%10s\n"%("btype", "k", "l0"))
topolfp.write("%10s%10.2f%10.3f\n"%("PP", 300.0, 0.65))
topolfp.close()

