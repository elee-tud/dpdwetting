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
    parser.add_option("-o", "--output", dest="output", default='moved.gro', type="str", help="Output file")
    parser.add_option("-v", "--velocity", dest="vel", default=0, type="float", help="New box size")
    (opts, args)=parser.parse_args()
    return opts.input, opts.output, opts.vel

inname, outname, dropvel=get_options()


wallgap=5.0
unitcell=0.5

title, topol, box, coord, vel=read_gro(inname)
for pvel in vel:
    pvel[2]+=dropvel
title="Droplet with v_z="+str(dropvel)
write_gro(outname, title, topol, box, coord, vel)
