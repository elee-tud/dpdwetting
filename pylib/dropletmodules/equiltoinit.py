#!/home/elee/programs/Python-3.6.10/bin/python

#------------------------------------------------------------------------------#
#  A program to put solvent molecules into a box randomly producing an input of
# a MDPD program.
#------------------------------------------------------------------------------#

import sys
from iocontrol.options import get_options as getopt
import subprocess





options=['-bx', '-by', '-bz', '-w', '-v', '-i', '-o', '-width', '-slope']
types=['float', 'float', 'float', 'float', 'float', 'str', 'str', 'float', 'float']
defaults=[80.0, 80, 60.0, 1.0, -2.0, 'final.gro', 'init.gro', 0, 0]
boxx, boxy, boxz, wallgap, velz, inname, outname, width, slope=getopt(sys.argv, options, types, defaults)

subprocess.call('centering_droplet.py -i '+inname, shell=True)
subprocess.call('move_droplet.py -i centered.gro -v '+str(velz), shell=True)
if slope==0.0:
    subprocess.call('add_wall.py -i movedrop.gro -x '+str(boxx)+' -y '+str(boxy)+' -z '+str(boxz)+' -o '+outname, shell=True)
else:
    subprocess.call('add_wiredwall.py -i movedrop.gro -x '+str(boxx)+' -y '+str(boxy)+' -z '+str(boxz)+' -o '+outname+' -width '+str(width)+' -slope '+str(slope), shell=True)


