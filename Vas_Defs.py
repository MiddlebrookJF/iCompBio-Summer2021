#!/usr/bin/python
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
    #this is the writhe of the entire chain takes the chain as input
    #uses gauss_lk
    wr=0.0
    acn=0.0
    nvertices=len(walk)
    if (nvertices<4):
        wr=0.0
        acn=0.0
    else:
        for j in range (1,nvertices-2,1):
            u1=np.array([walk[j-1]])
            u2=np.array([walk[j]])
            #u=u1-u2
            for i in range (j+2,nvertices,1):
                #print(i)
                #print(j)
                #print(walk[i-1])
                #print(walk[i])
                v1=np.array([walk[i-1]])
                v2=np.array([walk[i]])
                glk=gauss_lk(u1[0],u2[0],v1[0],v2[0])
                #print(u1[0],u2[0],v1[0],v2[0])
                #print(glk)
                wr=wr+2*glk
                acn=acn+2*abs(glk)
    return wr,acn

def randomvector(n):
    components = [np.random.normal() for i in range(n)]
    r = math.sqrt(sum(x*x for x in components))
    zv = [x/r for x in components]
    xv = np.random.randn(n)
    dot = xv.dot(zv)
    for r in range(0,len(zv)):
        xv[r] = xv[r] - (zv[r] * dot)
    xv /=  np.linalg.norm(xv)
    yv = np.cross(zv, xv)
    #print(np.array([xv, yv, zv]))
    #print(np.dot(xv,zv))
    #print(np.dot(xv, yv))
    #print(np.dot(yv, zv))
    return(np.array([xv,yv,zv]))


def vas_closed(verts):
    sum_vas=0
    all_vas=[]
    for p in range(0,500):
        #write as [[x,y,z],[x,y,z],...,[x,y,z]]
        v2=0
        nverts=len(verts)
        intersections=[]
        crossings=[]
        z=[]
        kcoords=[]
        ########################
        r=randomvector(3)
        for k in range(0,nverts):
            xk=((r[0][0]*verts[k][0])+(r[0][1]*verts[k][1])+(r[0][2]*verts[k][2]))/(math.sqrt((r[0][0]*r[0][0])+(r[0][1]*r[0][1])+(r[0][2]*r[0][2])))
            yk=((r[1][0]*verts[k][0])+(r[1][1]*verts[k][1])+(r[1][2]*verts[k][2]))/(math.sqrt((r[1][0]*r[1][0])+(r[1][1]*r[1][1])+(r[1][2]*r[1][2])))
            zk=((r[2][0]*verts[k][0])+(r[2][1]*verts[k][1])+(r[2][2]*verts[k][2]))/(math.sqrt((r[2][0]*r[2][0])+(r[2][1]*r[2][1])+(r[2][2]*r[2][2])))
            kcoord=[xk,yk,zk]
            kcoords.append(kcoord)
        #print(kcoords)
        ########################
        verts=kcoords ######### these are the projected coordinates
        if nverts<5:
            v2=0
        else:
            for j in range(0,nverts): 
                vect01=np.array(verts[j])
                if j<nverts-1:
                    vect02=np.array(verts[j+1])
                else:
                    vect02=np.array(verts[0])
                for i in range(j+2,nverts):
                    if i!=j+nverts-1:
                        vect11=np.array(verts[i])
                        if i<nverts-1:
                            vect12=np.array(verts[i+1])
                        else:
                            vect12=np.array(verts[0])
                        #print(j, vect01, vect02)
                        #print(i, vect11, vect12)
                        x01=vect01[0]
                        x02=vect02[0]
                        x11=vect11[0]
                        x12=vect12[0]
                        y01=vect01[1]
                        y02=vect02[1]
                        y11=vect11[1]
                        y12=vect12[1]
                        x0=x02-x01
                        x1=x12-x11
                        y0=y02-y01
                        y1=y12-y11
                        A=np.array([[x0,-x1],[y0,-y1]])
                        B=np.array([x11-x01,y11-y01])
                        if np.linalg.det(A)!=0:
                            X=np.linalg.inv(A).dot(B)
                            if (0<=X[0]<=1) and (0<=X[1]<=1):
                                I=np.array([X[0]*x0+x01,X[0]*y0+y01])
                                #print(I)
                                intersections.append(I)
                                zi0=vect01[2]+(X[0]*(vect02[2]-vect01[2]))
                                zi1=vect11[2]+(X[1]*(vect12[2]-vect11[2]))
                                #print(zi0, zi1)
                                crossings.append([j,i,X[0],X[1]])
                                z.append([zi0,zi1])
        #print(" ")
        #print(crossings)
        #print(" ")
        qcrossings=[]
        verts.append(verts[0])
        for k in range(0,len(crossings)-1):
            for l in range(k+1,len(crossings)):
                u=crossings[k][0]
                v=crossings[k][1]
                w=crossings[l][0]
                x=crossings[l][1]
                if (u<w) and (w<v) and (v<x) and ((z[k][0]>z[k][1] and z[l][0]<z[l][1]) or (z[k][0]<z[k][1] and z[l][0]>z[l][1])):
                    #print(crossings[k],crossings[l])
                    #print(u,v,w,x)
                    #print(verts[u])
                    #print(verts[u+1])
                    #print(verts[v])
                    #print(verts[v+1])
                    #print(np.sign(gauss_lk(np.array(verts[u]),np.array(verts[u+1]),np.array(verts[v]),np.array(verts[v+1]))))
                    v2=v2+float(np.sign(gauss_lk(np.array(verts[u]),np.array(verts[u+1]),np.array(verts[v]),np.array(verts[v+1]))))*float(np.sign(gauss_lk(np.array(verts[w]),np.array(verts[w+1]),np.array(verts[x]),np.array(verts[x+1]))))
                else:
                    if (u==w) and (crossings[k][2]<crossings[l][2]) and (w<v) and (v<x) and ((z[k][0]>z[k][1] and z[l][0]<z[l][1]) or (z[k][0]<z[k][1] and z[l][0]>z[l][1])):
                        #print("UW",crossings[k],crossings[l])
                        #print(u,v,w,x)
                        #print(verts[u])
                        #print(verts[u + 1])
                        #print(verts[v])
                        #print(verts[v + 1])
                        #print(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1]))))
                        v2 = v2 + float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]),np.array(verts[x + 1]))))
                    if (v==x) and (crossings[k][3]<crossings[l][3]) and (u<w) and (w<v) and ((z[k][0]>z[k][1] and z[l][0]<z[l][1]) or (z[k][0]<z[k][1] and z[l][0]>z[l][1])):
                        #print("VX",crossings[k],crossings[l])
                        #print(u,v,w,x)
                        #print(verts[u])
                        #print(verts[u + 1])
                        #print(verts[v])
                        #print(verts[v + 1])
                        #print(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1]))))
                        v2 = v2 + float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]),np.array(verts[x + 1]))))
                    if (v==w) and (crossings[l][2]<crossings[k][3]) and (u<w) and (v<x) and ((z[k][0]>z[k][1] and z[l][0]<z[l][1]) or (z[k][0]<z[k][1] and z[l][0]>z[l][1])):
                        #print("VW", crossings[k], crossings[l])
                        #print(u,v,w,x)
                        #print(verts[u])
                        #print(verts[u + 1])
                        #print(verts[v])
                        #print(verts[v + 1])
                        #print(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1]))))
                        #print(float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]), np.array(verts[x + 1])))))
                        v2 = v2 + float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]), np.array(verts[x + 1]))))
        #print("Second Vassiliev Invariant:")
        #print(v2)
        #print(" ")
        verts.pop(len(verts)-1)
        sum_vas+=v2
        all_vas.append(v2)
    print(all_vas)
    print(sum_vas)
    vas_avg=sum_vas/500
    print("Second Vassiliev:")
    print(vas_avg)

def vas_open(verts):
    if nverts<5:
            vas_avg=0
    else:
     sum_vas=0
     all_vas=[]
     for p in range(0,500):
        #write as [[x,y,z],[x,y,z],...,[x,y,z]]
        v2=0
        nverts=len(verts)
        intersections=[]
        crossings=[]
        z=[]
        kcoords=[]
        ########################
        r=randomvector(3)
        for k in range(0,nverts):
            xk=((r[0][0]*verts[k][0])+(r[0][1]*verts[k][1])+(r[0][2]*verts[k][2]))/(math.sqrt((r[0][0]*r[0][0])+(r[0][1]*r[0][1])+(r[0][2]*r[0][2])))
            yk=((r[1][0]*verts[k][0])+(r[1][1]*verts[k][1])+(r[1][2]*verts[k][2]))/(math.sqrt((r[1][0]*r[1][0])+(r[1][1]*r[1][1])+(r[1][2]*r[1][2])))
            zk=((r[2][0]*verts[k][0])+(r[2][1]*verts[k][1])+(r[2][2]*verts[k][2]))/(math.sqrt((r[2][0]*r[2][0])+(r[2][1]*r[2][1])+(r[2][2]*r[2][2])))
            kcoord=[xk,yk,zk]
            kcoords.append(kcoord)
        #print(kcoords)
        ########################
        verts=kcoords ######## these are the projected coordinates

        for j in range(0,nverts-4):
                vect01=np.array(verts[j])
                vect02=np.array(verts[j+1])
                for i in range(j+2,nverts-2):

                        vect11=np.array(verts[i])
                        vect12=np.array(verts[i+1])
                        #else:
                            #vect12=np.array(verts[0])
                        #print(j, vect01, vect02)
                        #print(i, vect11, vect12)
                        x01=vect01[0]
                        x02=vect02[0]
                        x11=vect11[0]
                        x12=vect12[0]
                        y01=vect01[1]
                        y02=vect02[1]
                        y11=vect11[1]
                        y12=vect12[1]
                        x0=x02-x01
                        x1=x12-x11
                        y0=y02-y01
                        y1=y12-y11
                        A=np.array([[x0,-x1],[y0,-y1]]) ###Kramer's rule for intersection 
                        B=np.array([x11-x01,y11-y01])
                        if np.linalg.det(A)!=0:
                            X=np.linalg.inv(A).dot(B)
                            if (0<=X[0]<=1) and (0<=X[1]<=1):
                                I=np.array([X[0]*x0+x01,X[0]*y0+y01])
                                #print(I)
                                intersections.append(I)
                                zi0=vect01[2]+(X[0]*(vect02[2]-vect01[2]))
                                zi1=vect11[2]+(X[1]*(vect12[2]-vect11[2]))
                                #print(zi0, zi1)
                                crossings.append([j,i,X[0],X[1]]) #### [edge j, edge I, 
                                z.append([zi0,zi1])
        #print(" ")
        #print(crossings)
        #print(" ")
        qcrossings=[]

        for k in range(0,len(crossings)-1):
            for l in range(k+1,len(crossings)):
                u=crossings[k][0]
                v=crossings[k][1]
                w=crossings[l][0]
                x=crossings[l][1]
                if (u<w) and (w<v) and (v<x) and ((z[k][0]>z[k][1] and z[l][0]<z[l][1]) or (z[k][0]<z[k][1] and z[l][0]>z[l][1])):
                    #print(crossings[k],crossings[l])
                    #print(u,v,w,x)
                    #print(verts[u])
                    #print(verts[u+1])
                    #print(verts[v])
                    #print(verts[v+1])
                    #print(np.sign(gauss_lk(np.array(verts[u]),np.array(verts[u+1]),np.array(verts[v]),np.array(verts[v+1]))))
                      v2=v2+float(np.sign(gauss_lk(np.array(verts[u]),np.array(verts[u+1]),np.array(verts[v]),np.array(verts[v+1]))))*float(np.sign(gauss_lk(np.array(verts[w]),np.array(verts[w+1]),np.array(verts[x]),np.array(verts[x+1]))))
                else:
                    if (u==w) and (crossings[k][2]<crossings[l][2]) and (w<v) and (v<x) and ((z[k][0]>z[k][1] and z[l][0]<z[l][1]) or (z[k][0]<z[k][1] and z[l][0]>z[l][1])):
                        #print("UW",crossings[k],crossings[l])
                        #print(u,v,w,x)
                        #print(verts[u])
                        #print(verts[u + 1])
                        #print(verts[v])
                        #print(verts[v + 1])
                        #print(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1]))))
                        v2 = v2 + float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]),np.array(verts[x + 1]))))
                    if (v==x) and (crossings[k][3]<crossings[l][3]) and (u<w) and (w<v) and ((z[k][0]>z[k][1] and z[l][0]<z[l][1]) or (z[k][0]<z[k][1] and z[l][0]>z[l][1])):
                        #print("VX",crossings[k],crossings[l])
                        #print(u,v,w,x)
                        #print(verts[u])
                        #print(verts[u + 1])
                        #print(verts[v])
                        #print(verts[v + 1])
                        #print(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1]))))
                        v2 = v2 + float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]),np.array(verts[x + 1]))))
                    if (v==w) and (crossings[l][2]<crossings[k][3]) and (u<w) and (v<x) and ((z[k][0]>z[k][1] and z[l][0]<z[l][1]) or (z[k][0]<z[k][1] and z[l][0]>z[l][1])):
                        #print("VW", crossings[k], crossings[l])
                        #print(u,v,w,x)
                        #print(verts[u])
                        #print(verts[u + 1])
                        #print(verts[v])
                        #print(verts[v + 1])
                        #print(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1]))))
                        #print(float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]), np.array(verts[x + 1])))))
                        v2 = v2 + float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]), np.array(verts[x + 1]))))
        #print("Second Vassiliev Invariant:")
        #print(v2)
        #print(" ")

        sum_vas+=v2
        all_vas.append(v2)
    #print(all_vas)
    #print(sum_vas)
    vas_avg=sum_vas/500
    #print("Second Vassiliev:")
    #print(vas_avg)
    print(vas_avg)



def vas_randwalk(n):
    #n = number of edges
    vks=[]
    vksum=0
    for t in range(0,500):
        walk=[[0,0,0]]
        #fig = plt.figure()
        #ax = plt.axes(projection='3d')
        for i in range (1,n+1):
            pt=[]
            r=randomvector(3)
            pt.append(walk[i-1][0]+r[2][0])
            pt.append(walk[i-1][1]+r[2][1])
            pt.append(walk[i-1][2]+r[2][2])
            walk.append(pt)
        print(walk)#
        vk=vas_open(walk)
        print(vk)#
        vks.append(vk)
        vksum+=vk
    vkavg=vksum/len(vks)
    print(vks)
    print(vkavg)

###########################

walk=[]
x=[]
y=[]
z=[]
for i in range (1,len(walk)+1):
    x.append(walk[i - 1][0])
    y.append(walk[i - 1][1])
    z.append(walk[i - 1][2])
fig=plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(x, y, z)
plt.show() 