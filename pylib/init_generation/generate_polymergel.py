#!/home/elee/.bin/python
##!/usr/bin/python

#------------------------------------------------------------------------------#
#  A program to put solvent molecules into a box randomly producing an input of
# a MDPD program.
#------------------------------------------------------------------------------#

import sys

sys.path.append("../")
import numpy as np
from iocontrol.options import get_options as getopt
import random as rd
from init_generation.coordinate import *
from init_generation.polymer_random_walk import *
from init_generation.pbc import *
from optparse import OptionParser

temp=1.0
mass=1.0
rd.seed()

parser=OptionParser()
parser.add_option("-n", "--total-number", dest="ntot", default=0, type="int", help="Total number of particles in a system")
parser.add_option("-d", "--polymer-density", dest="pd", default=0., type="float", help="Number density of polymer")
parser.add_option("-l", "--polymer-length", dest="pl", default=0., type="int", help="The length of a single polymer")
parser.add_option("-g", "--numpol-gel", dest="npg", default=0., type="int", help="The number of polymers in a gel")
parser.add_option("-c", "--crosslink-density", dest="ld", default=0., type="float", help="Crosslink density")
parser.add_option("-b", "--box-size", dest="box", default=0., type="float", help="Box size")
parser.add_option("-r", "--polymer-ratio", dest="pgbratio", default=0.2, type="float", help="The ratio of the box size that polymer located")
parser.add_option("-x", "--box-x", dest="boxx", default=0., type="float", help="Box size along x")
parser.add_option("-y", "--box-y", dest="boxy", default=0., type="float", help="Box size along y")
parser.add_option("-z", "--box-z", dest="boxz", default=0., type="float", help="Box size along z")
parser.add_option("-P", "--pol-pol", dest="pp", default=-20, type="float", help="Polymer-polymer attraction strength")
parser.add_option("-M", "--pol-sol", dest="ps", default=-30, type="float", help="Polymer-solvent attraction strength")
parser.add_option("-S", "--sol-sol", dest="ss", default=-20, type="float", help="Solvent-solvent atraction strength")
parser.add_option("-o", "--output", dest="outname", default="conf.gro", help="Output file name", metavar="FILE")
(opts, args)=parser.parse_args()

ntot, pd, pl, npg, ld, box, pgbratio, boxx, boxy, boxz, pp, ps, ss, outname=opts.ntot, opts.pd, opts.pl, opts.npg, opts.ld, opts.box, opts.pgbratio, opts.boxx, opts.boxy, opts.boxz, opts.pp, opts.ps, opts.ss, opts.outname
print(ntot)
"""
options=['-nt', '-pd', '-pl', '-ng', '-b', '-o' , '-r', '-x', '-y', '-z', '-pp', '-ss', '-ps'] 
types=['int', 'int', 'int', 'int', 'float', 'float', 'str', 'float', 'float', 'float', 'float', 'float', 'float', 'float']
defaults=[0, 0, 0, 1.0, 10.0, 'conf.gro', 0.2, 0., 0., 0., -40, -40, -40]
nsol, npol, lpol, bl, boxl, outname, polgenboxratio, boxx, boxy, boxz, pp, ss, ps=getopt(sys.argv, options, types, defaults)
if not boxx==0.:
    box=np.array([boxx, boxy, boxz])
else:
    box=np.array([boxl, boxl, boxl])
topol=[]
coord=[]
vel=[]
molidx=1
ptclidx=1
polgenbox=box*polgenboxratio
polgenmin=(box-polgenbox)/2

for i in range(npol):
    for j in range(lpol):
        vel.append(np.random.normal(0.0, temp/mass, 3))
        topol.append([molidx, 'POL', 'P', ptclidx])
        vel.append(np.random.normal(0.0, temp/mass, 3))
        ptclidx+=1
    else:
        polcoord=doublefolded(lpol, bl, box, polgenboxratio)
    coord+=polcoord

    molidx+=1
for i in range(len(coord)):
    coord[i]=get_particle_in_box(coord[i], box)
    
for i in range(nsol):
    topol.append([molidx, 'SOL', 'S', ptclidx])
    coord.append([rd.random()*box[0], rd.random()*box[1],rd.random()*box[2]])
    vel.append(np.random.normal(0.0, temp/mass, 3))
    ptclidx+=1
    molidx+=1

title="Polymer with solvent"
write_gro(outname, title, topol, box, coord, vel)
print(outname+" is generated.")

topolfp=open("topol.top", "w")
topolfp.write(";System information\n")
topolfp.write("[ System ]\n")
topolfp.write(";%9s%10s\n"%("Mol_name","Num_of_mols"))
if npol>0:
    topolfp.write("%10s%10d\n"%("POL", npol))
if nsol>0:
    topolfp.write("%10s%10d\n"%("SOL", nsol))

topolfp.write("\n")
topolfp.write(";Molecule Atom Information\n")
if nsol>0:
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
    if lpol>1:
        topolfp.write(";Molecule Bond Information\n")
        topolfp.write("[ POL Bonds ]\n")
        topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
        for i in range(lpol-1):
            topolfp.write("%10d%10d%10s\n"%(i, i+1, "PP"))
        if not isring==0:
            topolfp.write("%10d%10d%10s\n"%(lpol-1, 0, "PP"))

        topolfp.write("\n")
topolfp.write(";Nonbonded Parameters\n")
topolfp.write("[ Nonbonded ]\n")
topolfp.write(";%9s%10s%10s%10s%10s%10s\n"%("Atype1","Atype2","B", "r_B", "A", "r_A"))
if nsol>0:
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "S", 25, 0.75, ss, 1.0))
    if npol>0:
        topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "P", 25, 0.75, ps, 1.0))
        topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("P", "P", 25, 0.75, pp, 1.0))
topolfp.write("\n")
if lpol>1:
    topolfp.write(";Bond length Parameters\n")
    topolfp.write("[ Bondlength ]\n")
    topolfp.write(";%9s%10s%10s\n"%("btype", "k", "l0"))
    topolfp.write("%10s%10.2f%10.3f\n"%("PP", 300.0, 0.65))
topolfp.close()

"""
