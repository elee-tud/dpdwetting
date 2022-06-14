import random as rd
import numpy as np
import math
rd.seed()
def doublefolded(n, blength, box, startratio=1.0):
    coord=np.zeros((n, 3))
    halfn=int(n/2)
    if n%2==1:
        halfn+=1

    startbox=box*startratio
    center=np.zeros((halfn, 3))
    center[0]=np.array([rd.random()*startbox[0], rd.random()*startbox[1], rd.random()*startbox[2]])
    center[0]=center[0]+(box-startbox)/2.
#    print(center[0])

    for i in range(halfn-1):
        the=rd.random()*math.pi
        phi=rd.random()*math.pi*2.0
        dr=[blength*math.sin(the)*math.cos(phi), blength*math.sin(the)*math.sin(phi), blength*math.cos(the)]
        center[i+1]=center[i]+dr

    the=rd.random()*math.pi
    phi1=rd.random()*math.pi*2.0
    phi2=(phi1+math.pi)
    halfbl=blength/2

    for i in range(halfn):
        dr=[halfbl*math.sin(the)*math.cos(phi1), halfbl*math.sin(the)*math.sin(phi1), halfbl*math.cos(the)]
        coord[i]=center[i]+dr
        dr=[halfbl*math.sin(the)*math.cos(phi2), halfbl*math.sin(the)*math.sin(phi2), halfbl*math.cos(the)]
        coord[n-i-1]=center[i]+dr
    
    for i in range(n):
        if coord[i][0]>=box[0]:
            while coord[i][0]>=box[0]:
                coord[i][0]-=box[0]
        elif coord[i][0]<0:
            while coord[i][0]<0:
                coord[i][0]+=box[0]
        if coord[i][1]>=box[1]:
            while coord[i][1]>=box[1]:
                coord[i][1]-=box[1]
        elif coord[i][1]<0:
            while coord[i][1]<0:
                coord[i][1]+=box[1]
        if coord[i][2]>=box[2]:
            while coord[i][2]>=box[2]:
                coord[i][2]-=box[2]
        elif coord[i][2]<0:
            while coord[i][2]<0:
                coord[i][2]+=box[2]
        


    return coord.tolist()


