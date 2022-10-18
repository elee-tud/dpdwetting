#!/usr/bin/python

#------------------------------------------------------------------------------#
#  A program to put solvent molecules into a box randomly producing an input of
# a MDPD program.
#------------------------------------------------------------------------------#

import sys
sys.path.append('../')
import numpy as np
import random as rd
from coordinate import *
from pbc import *
from optparse import OptionParser

def get_options():
    parser=OptionParser()
    parser.add_option("-i", "--input", dest="input", default='nowall.gro', type="str", help="Input file")
    parser.add_option("-o", "--output", dest="output", default='conf.gro', type="str", help="Output file")
    parser.add_option("-x", "--box-size-x", dest="boxx", default=0, type="float", help="New box size along x")
    parser.add_option("-y", "--box-size-y", dest="boxy", default=0, type="float", help="New box size along y")
    parser.add_option("-z", "--box-size-z", dest="boxz", default=0, type="float", help="New box size along z")
    parser.add_option("-g", "--pillar-gap", dest="gap", default=0, type="float", help="Gap between pillars")
    parser.add_option("-t", "--height-gap", dest="height", default=0, type="float", help="Height of a pillar")
    parser.add_option("-w", "--width-gap", dest="width", default=0, type="float", help="Width of a pillar")
    parser.add_option("-d", "--pillar-direction", dest="direc", default=0, type="str", help="Direction of pillars")
    (opts, args)=parser.parse_args()
    return opts.input, opts.output, opts.boxx, opts.boxy, opts.boxz, opts.gap, opts.height, opts.width, opts.direc

inname, outname, boxx, boxy, boxz, gap, height, width, direc=get_options()



unitcell=0.5
wallgap=5.0

title, topol, box, coord, vel=read_gro(inname)
newbox=[boxx, boxy, boxz+2*wallgap+height]
if boxx!=0:
    newbox[0]=boxx
if boxy!=0:
    newbox[1]=boxy
if boxz!=0:
    newbox[2]=boxz+2*wallgap
diffx=(newbox[0]-box[0])/2
diffy=(newbox[1]-box[1])/2
diffz=(newbox[2]-box[2])/2
nboxx=int(newbox[0]/unitcell)
nboxy=int(newbox[1]/unitcell)
for i in range(len(coord)):
    coord[i][0]=get_particle_vector_in_box(coord[i][0], box[0])+diffx
    coord[i][1]=get_particle_vector_in_box(coord[i][1], box[1])+diffy
    coord[i][2]=get_particle_vector_in_box(coord[i][2], box[2])+diffz
last_molnum=topol[-1][0]
last_atnum=topol[-1][0]
index=1
ptclidx=len(coord)+1
molidx=topol[-1][0]+1


numwall=0
nsp=int((gap+width)/unitcell)
gapunit=int(gap/unitcell)
basenzh=3
botbasez=wallgap-(basenzh-1)*unitcell
topbasez=wallgap+boxz+height+(basenzh-1)*unitcell

#You have to write the code to add rough pillared walls to the droplet only system
for i in range(nboxx):
    for j in range(nboxy):
        if nsp==0:
            nzheight=basenzh
        elif i%nsp<gapunit or j%nsp<gapunit:
            nzheight=basenzh
        else:
            nzheight=int(height/unitcell)+basenzh
        
        for k in range(nzheight):
            coord.append([i*unitcell, j*unitcell, botbasez+k*unitcell])
            vel.append([0,0,0])
            topol.append([molidx, "WAL", "W", ptclidx])
            ptclidx+=1
            molidx+=1
            numwall+=1
#Surface at top
for i in range(nboxx):
    for j in range(nboxy):
        coord.append([i*unitcell, j*unitcell, topbasez])
        vel.append([0,0,0])
        topol.append([molidx, "WAL", "W", ptclidx])
        ptclidx+=1
        molidx+=1
        numwall+=1


title=title+" after rough wall added"
write_gro(outname, title, topol, newbox, coord, vel)
