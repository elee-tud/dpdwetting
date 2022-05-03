#!/home/elee/.bin/python
##!/p/software/juwels/stages/2022/software/Python/3.9.6-GCCcore-11.2.0/bin/python

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
options=['-s', '-p', '-l', '-bl',  '-o' , '-x', '-y', '-z', '-w', '-u', '-ix', '-iy', '-sw', '-gap', '-height', '-width', '-dir']
types=['int', 'int', 'int', 'float', 'str', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'float', 'str']
defaults=[0, 0, 0, 0.65, 'conf.gro',  0., 0., 0., 5.0, 0.5, 0., 0., -10, 0, 0, 0, 'x']
nsol, npol, lpol, bl, outname, boxx, boxy, boxz, wallgap, unitcell, inboxx, inboxy, solv_wall, gap, height, width, direc=getopt(sys.argv, options, types, defaults)


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
filmz=boxz-height-2*unitcell
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
        polcoord[j][2]+=wallgap+height+unitcell

    coord+=polcoord
    molidx+=1

for i in range(nsol):
    topol.append([molidx, 'SOL', 'S', ptclidx])
    coord.append([rd.random()*inboxx+(boxx-inboxx)/2, rd.random()*inboxy+(boxy-inboxy)/2,rd.random()*filmz+wallgap+height+unitcell])
    vel.append(np.random.normal(0.0, temp/mass, 3))
    ptclidx+=1
    molidx+=1

numwall=0
nsp=int((gap+width)/unitcell)
gapunit=int(gap/unitcell)
basenzh=3
botbasez=wallgap-(basenzh-1)*unitcell
topbasez=wallgap+boxz+height+(basenzh-1)*unitcell

if direc=='x':
    for i in range(nboxx):
        if nsp==0:
            nzheight=basenzh
        elif i%nsp<gapunit:
            nzheight=basenzh
        else:
            nzheight=int(height/unitcell)+basenzh
        for j in range(nboxy):
            for k in range(nzheight):
                coord.append([i*unitcell, j*unitcell, botbasez+k*unitcell])
                vel.append([0,0,0])
                topol.append([molidx, "WAL", "W", ptclidx])
                ptclidx+=1
                molidx+=1
                numwall+=1
    for i in range(nboxx):
        if nsp==0:
            nzheight=basenzh
        elif i%nsp<gapunit:
            nzheight=basenzh
        else:
            nzheight=int(height/unitcell)+basenzh
        for j in range(nboxy):
            for k in range(nzheight):
                coord.append([i*unitcell, j*unitcell, topbasez-k*unitcell])
                vel.append([0,0,0])
                topol.append([molidx, "RAL", "W", ptclidx])
                ptclidx+=1
                molidx+=1
elif direc=='y':
    for i in range(nboxy):
        if i%nsp<gapunit:
            nzheight=basenzh
        else:
            nzheight=int(height/unitcell)+basenzh
        for j in range(nboxx):
            for k in range(nzheight):
                coord.append([j*unitcell, i*unitcell, botbasez+k*unitcell])
                vel.append([0,0,0])
                topol.append([molidx, "WAL", "W", ptclidx])
                ptclidx+=1
                molidx+=1
                numwall+=1
    for i in range(nboxy):
        if i%nsp<gapunit:
            nzheight=basenzh
        else:
            nzheight=int(height/unitcell)+basenzh
        for j in range(nboxx):
            for k in range(nzheight):
                coord.append([j*unitcell, i*unitcell, topbasez-k*unitcell])
                vel.append([0,0,0])
                topol.append([molidx, "RAL", "W", ptclidx])
                ptclidx+=1
                molidx+=1

            

box[2]=box[2]+2*wallgap+height

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
topolfp.write("%10s%10d%10s\n"%("WAL", numwall , "frozen"))
topolfp.write("%10s%10d%10s\n"%("RAL", numwall, "frozen"))

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

