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
from dropletmodules.polymer_random_walk import *
from dropletmodules.pbc import *

temp=0.5
mass=1.0
rd.seed()
box=[0.0, 0.0, 0.0]
options=['-s', '-p', '-l', '-w', '-bl', '-b', '-o', '-pw', '-pp', '-ps', '-ss']
types=['int', 'int', 'int', 'int', 'float', 'float', 'str', 'float', 'float', 'float', 'float']
defaults=[0, 0, 0, 0, 0.65, 80.0, 'topol.top', -30, -40, -40, -40]
nsol, npol, lpol, nwal, bl, boxx, outname, pwparm, pp, ps, ss=getopt(sys.argv, options, types, defaults)
if nwal==0:
    nwal=2*(boxx/0.5)*(boxx/0.5)
topolfp=open(outname, "w")
topolfp.write(";System information\n")
topolfp.write("[ System ]\n")
topolfp.write(";%9s%10s\n"%("Mol_name","Num_of_mols"))
if npol>0:
    topolfp.write("%10s%10d\n"%("POL", npol))
if nsol>0:
    topolfp.write("%10s%10d\n"%("SOL", nsol))
if nwal>0:
    topolfp.write("%10s%10d%10s\n"%("WAL", nwal, "frozen"))

topolfp.write("\n")
topolfp.write(";Molecule Atom Information\n")
if nsol>0:
    topolfp.write("[ SOL Atoms ]\n")
    topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
    topolfp.write("%10d%10s%10.1f\n"%(0, "W", 1.0))
    topolfp.write("\n")
if nwal>0:
    topolfp.write("[ WAL Atoms ]\n")
    topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
    topolfp.write("%10d%10s%10.1f\n"%(0, "S", 1.0))
    topolfp.write("\n")
if npol>0:
    topolfp.write("[ POL Atoms ]\n")
    topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
    for i in range(lpol):
        topolfp.write("%10d%10s%10.1f\n"%(i, "P", 1.0))
    topolfp.write("\n")
if lpol>1:
    topolfp.write(";Molecule Bond Information\n")
    topolfp.write("[ POL Bonds ]\n")
    topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
    for i in range(lpol-1):
        topolfp.write("%10d%10d%10s\n"%(i, i+1, "PP"))
    topolfp.write("\n")
topolfp.write(";Nonbonded Parameters\n")
topolfp.write("[ Nonbonded ]\n")
topolfp.write(";%9s%10s%10s%10s%10s%10s\n"%("Atype1","Atype2","B", "r_B", "A", "r_A"))
if nsol>0:
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("W", "W", 25, 0.75, ss, 1.0))
    if npol>0:
        topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("W", "P", 25, 0.75, ps, 1.0))
if nwal>0:
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("W", "S", 25, 0.75, -10, 1.0))
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "S", 25, 0.75, -40, 1.0))
    if npol>0:
        topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("P", "S", 25, 0.75, pwparm, 1.0))
if npol>0:
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("P", "P", 25, 0.75, pp, 1.0))
topolfp.write("\n")
if lpol>1:
    topolfp.write(";Bond length Parameters\n")
    topolfp.write("[ Bondlength ]\n")
    topolfp.write(";%9s%10s%10s\n"%("btype", "k", "l0"))
    topolfp.write("%10s%10.2f%10.3f\n"%("PP", 300.0, 0.65))
topolfp.close()


