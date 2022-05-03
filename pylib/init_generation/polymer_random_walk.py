import random as rd
import numpy as np
import math
rd.seed()
def rdpolymer(n, blength, box, first=False, boundary=[False, False, False]):
    coord=np.zeros((n, 3))
    if not first:
        coord[0]=np.array([rd.random()*box[0], rd.random()*box[1], rd.random()*box[2]])
    else:
        coord[0]=first
    for i in range(n-1):
        the=rd.random()*math.pi
        phi=rd.random()*math.pi*2.0
        dr=[blength*math.sin(the)*math.cos(phi), blength*math.sin(the)*math.sin(phi), blength*math.cos(the)]
        r=coord[i]+dr
        if boundary[0]:
            if r[0]<0:
                r[0]=-r[0]
            elif r[0]>=box[0]:
                r[0]=2*box[0]-r[0]
        if boundary[1]:
            if r[1]<0:
                r[1]=-r[1]
            elif r[1]>=box[1]:
                r[1]=2*box[1]-r[1]
        if boundary[2]:
            if r[2]<0:
                r[2]=-r[2]
            elif r[2]>=box[2]:
                r[2]=2*box[2]-r[2]
        coord[i+1]=r


    return coord.tolist()

def branchedpolymer(lbb, dbr, lbr, blength, box):
    nbr=int(lbb/dbr)
    ntot=lbb+nbr*lbr
    coord=np.zeros((ntot, 3))
    coord[0]=np.array([rd.random()*box[0], rd.random()*box[1], rd.random()*box[2]])

    for i in range(lbb-1):
        the=rd.random()*math.pi
        phi=rd.random()*math.pi*2.0
        dr=[blength*math.sin(the)*math.cos(phi), blength*math.sin(the)*math.sin(phi), blength*math.cos(the)]
        coord[i+1]=coord[i]+dr

    for i in range(nbr):
        bpoint=(i+1)*dbr-1
        the=rd.random()*math.pi
        phi=rd.random()*math.pi*2.0
        dr=[blength*math.sin(the)*math.cos(phi), blength*math.sin(the)*math.sin(phi), blength*math.cos(the)]
        coord[lbb+i*lbr]=coord[bpoint]+dr
        for j in range(lbr-1):
            the=rd.random()*math.pi
            phi=rd.random()*math.pi*2.0
            dr=[blength*math.sin(the)*math.cos(phi), blength*math.sin(the)*math.sin(phi), blength*math.cos(the)]
            coord[lbb+i*lbr+j+1]=coord[lbb+i*lbr+j]+dr

    return coord.tolist()
            





