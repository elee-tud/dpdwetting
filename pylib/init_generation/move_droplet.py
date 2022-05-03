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
from dropletmodules.pbc import *



wallgap=5.0
unitcell=0.5
options=['-v', '-i', '-o']
types=['float', 'str', 'str']
defaults=[0.0, 'nowall.gro', 'movedrop.gro']
dropvel, inname, outname=getopt(sys.argv, options, types, defaults)

title, topol, box, coord, vel=read_gro(inname)
for pvel in vel:
    pvel[2]+=dropvel
title="Droplet with v_z="+str(dropvel)
write_gro(outname, title, topol, box, coord, vel)
