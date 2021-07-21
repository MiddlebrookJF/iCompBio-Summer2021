import os
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#mpl.use('Agg')
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D

def gauss_lk(a1,a2,b1,b2): #this is the Gauss linking number of two vectors
    #print("the vectors are")
    #print(a1)
    #print(a2)
    #print(b1)
    #print(b2)
    ra=a2-a1
    rb=b2-b1
    r00=a1-b1
    r01=a1-b2
    r10=a2-b1
    r11=a2-b2
    v1=np.cross(r00,r01)
    v1=v1/np.linalg.norm(v1)
    v2=np.cross(r01,r11)
    v2=v2/np.linalg.norm(v2)
    v3=np.cross(r11,r10)
    v3=v3/np.linalg.norm(v3)
    v4=np.cross(r10,r00)
    v4=v4/np.linalg.norm(v4)
    d1=np.dot(v1,v2)
    d2=np.dot(v2,v3)
    d3=np.dot(v3,v4)
    d4=np.dot(v4,v1)
    as1=np.arcsin(d1)
    as2=np.arcsin(d2)
    as3=np.arcsin(d3)
    as4=np.arcsin(d4)
    aux1=np.cross(ra,rb)
    aux=np.dot(aux1,r00)
    alk=np.sign(aux)*(as1+as2+as3+as4)/(4*math.pi)
    #sgn=np.sign(aux)

    return alk

def compute_wr(walk):
    #this is the writhe of the entire chain, takes the chain as input
    #uses gauss_lk
    wr=0.0
    acn=0.0
    nvertices=len(walk)
    if (nvertices<4):
        wr=0.0
        acn=0.0
    else:
        for j in range (1, nvertices-2, 1):
            u1 = np.array( [walk[j-1]] )
            u2 = np.array( [walk[j]] )
            #u=u1-u2
            for i in range (j+2,nvertices,1):
                #print(i)
                #print(j)
                #print(walk[i-1])
                #print(walk[i])
                v1 = np.array( [walk[i-1]] )
                v2 = np.array( [walk[i]] )
                glk = gauss_lk( u1[0], u2[0], v1[0], v2[0] )
                #print(u1[0],u2[0],v1[0],v2[0])
                #print(glk)
                wr = wr + 2 *glk
                acn = acn + 2 * abs(glk)
    return wr,acn
