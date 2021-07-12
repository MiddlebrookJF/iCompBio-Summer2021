#!/usr/bin/python
import os
#import MDAnalysis
#import MDAnalysis.coordinates
import matplotlib as mpl
mpl.use('Agg')
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D
import math
import numpy as np
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
import pickle as pl
from operator import itemgetter
import scipy.stats

from matplotlib.colors import BoundaryNorm


def gauss_lk(a1,a2,b1,b2): #this is the Gauss linking number of two vectors
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
    return alk



def compute_lk(walk1,walk2):
    #print(walk1)
    #print(walk2)
    #this is the linking number between the walk (x,y,z) and the walk (x1,y1,z1) translated by (m1,m2,m3)
    #uses gauss_lk
    lk=0
    icn=0
    for j in range (1,len(walk1),1):
        u1=np.array([walk1[j-1]])
        u2=np.array([walk1[j]])
        #u=u1-u2
        for i in range (1,len(walk2),1):
                v1=np.array([walk2[i-1]])
                v2=np.array([walk2[i]])
                glk=gauss_lk(u1[0],u2[0],v1[0],v2[0])
                lk=lk+glk
                icn=icn+abs(glk)
    return lk    #, icn



def linking_fmatrix(walk):
    walk1=[]
    walk2=[]
    walk3=[]
    lkm=[[0 for k in range(len(walk))] for j in range(len(walk))]
    #lkm=np.zeros((len(walk),len(walk)))
    walk1.append(walk[0])
    for i in range(2,len(walk)-4,1):
        #print("i")
        #print(i)
        walk1.append(walk[i-1])
        walk2=[]
        for j in range (i+1,len(walk)-2,1):
            #print("j")
            #print(j)
            walk3=[]
            walk2.append(walk[j-1])
            if len(walk2)>2:
               for k in range (j+1,len(walk)+1,1):
                   walk3.append(walk[k-1])
               #print("WALK1")
               #print(walk1)
               #print("WALK2")
               #print(walk2)
               lkm[i-1][j-1]=compute_lk(walk1,walk2)
               #print(lkm[i-1][j-1])
               #print("WALK3")
               #print(walk3)
               lkm[j-1][i-1]=compute_lk(walk2,walk3)
               #print(lkm[j-1][i-1])
    return lkm



def lk_fprint(lkfv1,protl): 
    print("H00")
    fig=plt.figure()
    print("H001")
    ax=fig.add_subplot(111)
    print("H0")
    ax.set_aspect('equal')
    cmap = plt.get_cmap('PuOr')
    # define the colormap
    cmap = plt.get_cmap('PuOr')
    # extract all colors from the .jet map
    cmaplist = [cmap(i) for i in range(cmap.N)]
    # create the new map
    cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N)
    print("This is it")
    print(lkfv1)
    lkfv=np.array(lkfv1)
    print("The min is")
    print(np.min(lkfv))
    lkfv2=np.matrix(lkfv)
    print(lkfv2)
    #print("H00")
    #fig=plt.figure()
    #print("H001")
    #ax=fig.add_subplot(111)
    #print("H0")
    #ax.set_aspect('equal')
    vm=max(abs(np.min(lkfv)),abs(np.max(lkfv)))
    norm=mpl.colors.Normalize(vmin=-vm,vmax=vm,clip=False)
    print("H1")
    plt.imshow(lkfv2,interpolation='none',cmap=cmap,norm=norm) #cmap=plt.cm.ocean)
    print("H2")
    plt.colorbar()
    #plt.show()
    plt.savefig("draft"+".png")
    #min_val, max_val = 0, 1
    #for i in xrange(protl):
    #    for j in xrange(protl):
    #        c = lkfv[i][j]
    #        ax.text(i+0.5, j+0.5, str(c), va='center', ha='center')
    #plt.matshow(lkfv, cmap=plt.cm.Blues)
    #ax.set_xlim(min_val, max_val)
    #ax.set_ylim(min_val, max_val)
    #ax.set_xticks(np.arange(max_val))
    #ax.set_yticks(np.arange(max_val))
    #ax.grid()



def read_coord(prot1,start,end): #,start,end): #period,tk,ampl
#this will read the file for a given period and then print it in the form to be read by lkp code
    prot=prot1 # +".pdb.n"; #"_dt2.pdb";
    print("This is protein")
    print(prot1)
    #h1=prot1+"_wr"
    f = open(prot,'r');
    l=0
    walk=[]
    k=0
    for line in f:
        k=k+1
        b=line.split()
        print(b)
        if (k>start-1) and (k<end+1):
           walk.append([float(b[0]),float(b[1]),float(b[2])])
    protl=end-start
    print("The length of the protein is:")
    print(protl)
    return walk,protl


def read_all_coord():
    for i in range (1,2,1): #len(protlist)+1,1): #len(protlist)/3
        (walk,protl)=read_coord("2fma_xyz2",1,59) #   str(protlist[i-1]),start,end) 
        lkfv=linking_fmatrix(walk)
        lk_fprint(lkfv,protl)


read_all_coord()
