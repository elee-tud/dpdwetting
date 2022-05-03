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
options=['-s', '-p', '-l', '-bl',  '-o' , '-x', '-y', '-z', '-w', '-u', '-ix', '-iy', '-sw']
types=['int', 'int', 'int', 'float', 'str', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float']
defaults=[0, 0, 0, 0.65, 'conf.gro',  0., 0., 0., 5.0, 0.5, 0., 0., -10]
nsol, npol, lpol, bl, outname, boxx, boxy, boxz, wallgap, unitcell, inboxx, inboxy, solv_wall=getopt(sys.argv, options, types, defaults)


box=np.array([boxx, boxy, boxz])
if inboxx==0.:
    inboxx=boxx
if inboxy==0.:
    inboxy=boxy

topol=[]
coord=[]
vel=[]
molidx=1
ptclidx=1
nboxx=int(box[0]/unitcell)
nboxy=int(box[1]/unitcell)
filmz=boxz-2*unitcell
boxfilm=np.array([inboxx, inboxy, filmz])

for i in range(npol):
    for j in range(lpol):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'POL', 'P', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    if inboxx==0. and inboxy==0.:
        polcoord=rdpolymer(lpol, bl, boxfilm, boundary=[False, False, True])
    elif inboxx==boxx and inboxy!=boxy:
        polcoord=rdpolymer(lpol, bl, boxfilm, boundary=[False, True, True])
    elif inboxx!=boxx and inboxy==boxy:
        polcoord=rdpolymer(lpol, bl, boxfilm, boundary=[True, False, True])
    elif inboxx!=boxx and inboxy!=boxy:
        polcoord=rdpolymer(lpol, bl, boxfilm, boundary=[True, True, True])
    for j in range(lpol):
        polcoord[j][0]+=(boxx-inboxx)/2
        polcoord[j][1]+=(boxy-inboxy)/2
        polcoord[j][2]+=wallgap+unitcell

    coord+=polcoord
    molidx+=1

for i in range(nsol):
    topol.append([molidx, 'SOL', 'W', ptclidx])
    coord.append([rd.random()*inboxx+(boxx-inboxx)/2, rd.random()*inboxy+(boxy-inboxy)/2,rd.random()*filmz+wallgap+unitcell])
    vel.append(np.random.normal(0.0, temp/mass, 3))
    ptclidx+=1
    molidx+=1

for i in range(nboxx):
    for j in range(nboxy):
        coord.append([i*unitcell, j*unitcell, wallgap])
        vel.append([0,0,0])
        topol.append([molidx, "WAL", "S", ptclidx])
        ptclidx+=1
        molidx+=1
for i in range(nboxx):
    for j in range(nboxy):
        coord.append([i*unitcell, j*unitcell, boxz+wallgap])
        vel.append([0,0,0])
        topol.append([molidx, "RAL", "R", ptclidx])
        ptclidx+=1
        molidx+=1

            

box[2]=box[2]+2*wallgap 

title="polymer grafting wall"
write_gro(outname, title, topol, box, coord, vel)
print(outname+" is generated.")

topolfp=open("topol.top", "w")
topolfp.write(";System information\n")
topolfp.write("[ System ]\n")
topolfp.write(";%9s%10s\n"%("Mol_name","Num_of_mols"))
if npol>0:
    topolfp.write("%10s%10d\n"%("POL", npol))
topolfp.write("%10s%10d\n"%("SOL", nsol))
topolfp.write("%10s%10d%10s\n"%("WAL", nboxx*nboxy , "frozen"))
topolfp.write("%10s%10d%10s\n"%("RAL", nboxx*nboxy, "frozen"))

topolfp.write("\n")
topolfp.write(";Molecule Atom Information\n")
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
    topolfp.write(";Molecule Bond Information\n")
    topolfp.write("[ POL Bonds ]\n")
    topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
    for i in range(lpol-1):
        topolfp.write("%10d%10d%10s\n"%(i, i+1, "PP"))
    topolfp.write("\n")

topolfp.write("[ WAL Atoms ]\n")
topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
topolfp.write("%10d%10s%10.1f\n"%(0, "W", 1.0))
topolfp.write("\n")
topolfp.write("[ RAL Atoms ]\n")
topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
topolfp.write("%10d%10s%10.1f\n"%(0, "W", 1.0))
topolfp.write("\n")

topolfp.write(";Nonbonded Parameters\n")
topolfp.write("[ Nonbonded ]\n")
topolfp.write(";%9s%10s%10s%10s%10s%10s\n"%("Atype1","Atype2","B", "r_B", "A", "r_A"))
topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "S", 25, 0.75, -40, 1.0))
topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "W", 25, 0.75, solv_wall, 1.0))
topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("W", "W", 25, 0.75, -40, 1.0))
if npol>0:
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "P", 25, 0.75, -40, 1.0))
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("P", "P", 25, 0.75, -40, 1.0))
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("P", "W", 25, 0.75, -30, 1.0))
topolfp.write("\n")

if lpol>1:
    topolfp.write(";Bond length Parameters\n")
    topolfp.write("[ Bondlength ]\n")
    topolfp.write(";%9s%10s%10s\n"%("btype", "k", "l0"))
    topolfp.write("%10s%10.2f%10.3f\n"%("PP", 300.0, 0.65))
topolfp.close()

