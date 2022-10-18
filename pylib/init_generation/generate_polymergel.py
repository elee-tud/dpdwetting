#!/home/elee/.bin/python
##!/usr/bin/python

#------------------------------------------------------------------------------#
#  A program to put solvent molecules into a box randomly producing an input of
# a MDPD program.
#------------------------------------------------------------------------------#

import sys
sys.path.append('../')

import numpy as np
import random as rd
from init_generation.coordinate import write_gro
from init_generation.polymer_random_walk import rdpolymer
# from init_generation.pbc import *
from optparse import OptionParser
import math

temp=1.0
mass=1.0
rd.seed()

def get_options():
    parser=OptionParser()
    parser.add_option("-n", "--total-number", dest="ntot", default=0, type="int", help="Total number of particles in a system")
    parser.add_option("-d", "--polymer-density", dest="pd", default=0., type="float", help="Number density of polymer")
    parser.add_option("-l", "--polymer-length", dest="pl", default=0., type="int", help="The length of a single polymer")
    parser.add_option("-g", "--numpol-gel", dest="npg", default=0., type="int", help="The number of polymers in a gel")
    parser.add_option("-c", "--crosslink-density", dest="ld", default=0., type="float", help="Crosslink density")
    parser.add_option("-b", "--box-size", dest="boxl", default=0., type="float", help="Box size")
    parser.add_option("-s", "--lattice-size", dest="lsize", default=0., type="float", help="Lattice size of the gel placed")
    parser.add_option("-t", "--distance-neighbor", dest="ndc", default=2.0, type="float", help="Distance criterion for the neighbor")
    parser.add_option("-P", "--pol-pol", dest="pp", default=-20, type="float", help="Polymer-polymer attraction strength")
    parser.add_option("-M", "--pol-sol", dest="ps", default=-30, type="float", help="Polymer-solvent attraction strength")
    parser.add_option("-S", "--sol-sol", dest="ss", default=-40, type="float", help="Solvent-solvent atraction strength")
    parser.add_option("-o", "--output", dest="outname", default="conf.gro", help="Output file name", metavar="FILE")
    (opts, args)=parser.parse_args()
    return opts.ntot, opts.pd, opts.pl, opts.npg, opts.ld, opts.boxl, opts.pp, opts.ps, opts.ss, opts.outname, opts.lsize, opts.ndc

def neighbor_list(coordinate, len_pol, distcrit=1.5):
    coordinate=np.array(coordinate)
    n_ptcls=len(coordinate)
    nblist=[[] for i in range(n_ptcls)]
    for i in range(n_ptcls):
        for j in range(i+1,n_ptcls):
            if np.linalg.norm(coordinate[i]-coordinate[j])<distcrit and int(i/len_pol) != int(j/len_pol):
                nblist[i].append(j)
                nblist[j].append(i)
            
    return nblist
        
        
    

n_tot_ptcls, pol_dens, pol_len, n_pol_gel, cl_dens, boxl, att_pp, att_ps, att_ss, outname, lat_gel_size, distcrit=get_options()

#box size#
box=np.array([boxl, boxl, boxl])
topol=[]
coord=[]
vel=[]
molidx=1
ptclidx=1
n_pol_ptcls=int(n_tot_ptcls*pol_dens)
n_gel=int(n_pol_ptcls/pol_len/n_pol_gel)
n_pol_ptcls_gel=int(n_pol_ptcls/n_gel)
n_cl_per_pol=int(cl_dens*pol_len)
n_pol_gel=int(n_pol_ptcls_gel/pol_len)
if lat_gel_size==0.: lat_gel_size=n_pol_ptcls_gel**0.5*0.5

lat_num=math.ceil(n_gel**(1/3))
lat_box_size=boxl/lat_num
if lat_box_size<lat_gel_size:
    raise ValueError("Lattice size of the gel to be placed is larger than the lattice size of the box determined by the number of gels.")
n_sol=int(n_tot_ptcls*(1-pol_dens))

for i in range(n_gel):
    xidx=int(i/(lat_num**2))
    yidx=int((i%(lat_num**2))/lat_num)
    zidx=(i%(lat_num**2))%lat_num
    genbox_origin=[(xidx+0.5)*lat_box_size-lat_gel_size*0.5, (yidx+0.5)*lat_box_size-lat_gel_size*0.5, (zidx+0.5)*lat_box_size-lat_gel_size*0.5]
    gelcoord=[]
    for j in range(n_pol_gel): 
        gelcoord+=rdpolymer(pol_len, 0.65, [lat_gel_size, lat_gel_size, lat_gel_size], first=False, boundary=[True, True, True])
        topol+=[[i+1, 'POL', "P", i*(n_pol_gel+pol_len)+j*pol_len+k+1] for k in range(pol_len)]
    nblist=neighbor_list(gelcoord, pol_len, distcrit)
    n_cl=0
    crosslinks=[]
    for i in range(n_pol_gel):
        for j in range(n_cl_per_pol):
            trynum=0
            while trynum < 100:
                pick=int(rd.random()*pol_len)+pol_len*i
                if len(nblist[pick])>0:
                    nbpick=int(rd.random()*len(nblist[pick]))
                    crosslinks.append([pick, nblist[pick][nbpick]])
                    break
                trynum+=1
            if trynum==100: raise ValueError("Too many failures to create cross links.")
    for c in gelcoord:
        coord.append([c[0]+genbox_origin[0], c[1]+genbox_origin[1], c[2]+genbox_origin[2]])
    molidx+=1
    ptclidx+=n_pol_ptcls_gel            
                
        
            
            
        
    
 
vel+=[np.random.normal(0.0, temp/mass, 3) for i in range(n_pol_ptcls)]


    

    
for i in range(n_sol):
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
if n_gel>0:
    topolfp.write("%10s%10d\n"%("POL", n_gel))
if n_sol>0:
    topolfp.write("%10s%10d\n"%("SOL", n_sol))

topolfp.write("\n")
topolfp.write(";Molecule Atom Information\n")
if n_sol>0:
    topolfp.write("[ SOL Atoms ]\n")
    topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
    topolfp.write("%10d%10s%10.1f\n"%(0, "S", 1.0))
    topolfp.write("\n")
if n_gel>0:
    topolfp.write("[ POL Atoms ]\n")
    topolfp.write(";%9s%10s%10s\n"%("Aindex", "Atype", "mass"))
    for i in range(n_pol_ptcls_gel):
        topolfp.write("%10d%10s%10.1f\n"%(i, "P", 1.0))
    topolfp.write("\n")
    if pol_len>1:
        topolfp.write(";Molecule Bond Information\n")
        topolfp.write("[ POL Bonds ]\n")
        topolfp.write(";%9s%10s%10s\n"%("index1", "index2", "btype"))
        for i in range(n_pol_gel):
            for j in range(pol_len-1):
                topolfp.write("%10d%10d%10s\n"%(i*pol_len+j, i*pol_len+j+1, "PP"))
        for cl in crosslinks:
            topolfp.write("%10d%10d%10s\n"%(cl[0], cl[1], "PP"))
            

        topolfp.write("\n")
        
topolfp.write(";Nonbonded Parameters\n")
topolfp.write("[ Nonbonded ]\n")
topolfp.write(";%9s%10s%10s%10s%10s%10s\n"%("Atype1","Atype2","B", "r_B", "A", "r_A"))
if n_sol>0:
    topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "S", 25, 0.75, att_ss, 1.0))
    if n_gel>0:
        topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("S", "P", 25, 0.75, att_ps, 1.0))
        topolfp.write("%10s%10s%10.2f%10.2f%10.2f%10.2f\n"%("P", "P", 25, 0.75, att_pp, 1.0))
topolfp.write("\n")
if pol_len>1:
    topolfp.write(";Bond length Parameters\n")
    topolfp.write("[ Bondlength ]\n")
    topolfp.write(";%9s%10s%10s\n"%("btype", "k", "l0"))
    topolfp.write("%10s%10.2f%10.3f\n"%("PP", 300.0, 0.65))
topolfp.close()


