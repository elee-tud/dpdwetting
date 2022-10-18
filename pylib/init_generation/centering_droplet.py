#!/usr/bin/python

#------------------------------------------------------------------------------#
#  A program to put solvent molecules into a box randomly producing an input of
# a MDPD program.
#------------------------------------------------------------------------------#

import sys
sys.path.append('../')
import numpy as np
from optparse import OptionParser
import random as rd
from coordinate import *
from pbc import *



def get_options():
    parser=OptionParser()
    parser.add_option("-i", "--input", dest="input", default='before.gro', type="str", help="Input file")
    parser.add_option("-o", "--output", dest="output", default='centered.gro', type="str", help="Output file")
    parser.add_option("-b", "--box size", dest="boxl", default=-1, type="float", help="New box size")
    (opts, args)=parser.parse_args()
    return opts.input, opts.output, opts.boxl

inname, outname, nbox=get_options()

title, topol, box, coord, vel=read_gro(inname)
if nbox==-1.0:
    newbox=box
else:
    newbox=np.array([nbox, nbox, nbox])

refrange=[[box[0]/2-1.0, box[0]/2+1.0],[box[1]/2-1.0, box[1]/2+1.0],[box[2]/2-1.0, box[2]/2+1.0]]
ref=np.array(coord[1])
for i in range(len(coord)):
    if coord[i][0]>=refrange[0][0] and coord[i][0]<refrange[0][1] and coord[i][1]>=refrange[1][0] and coord[i][1]<refrange[1][1] and coord[i][2]>=refrange[2][0] and coord[i][2]<refrange[2][1]:
        ref=np.array(coord[i])
        break

com=np.array([0.0,0.0,0.0])
for i in range(len(coord)):
    rij=get_shortest_image_vector(np.array(coord[i])-ref, box)
    com=com+rij


    

com=com/len(coord)+ref

for i in range(len(coord)):
    coord[i]=get_particle_in_box(np.array(coord[i])-com+np.array(newbox)/2, np.array(newbox))

title=title+" after centered"
write_gro(outname, title, topol, newbox, coord, vel)
