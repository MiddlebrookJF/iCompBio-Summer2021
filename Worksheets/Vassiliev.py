import os
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
from multiprocessing import Pool

#Gauss linking integral of two vectors given their endpoints, all inputs are np.arrays
def gauss_lk(a1,a2,b1,b2): 
    #calculating the sides of the quadrilateral
    ra=a2-a1
    rb=b2-b1
    r00=a1-b1
    r01=a1-b2
    r10=a2-b1
    r11=a2-b2

    #calculate and normalize side cross products
    v1=np.cross(r00,r01)
    v1=v1/np.linalg.norm(v1)
    v2=np.cross(r01,r11)
    v2=v2/np.linalg.norm(v2)
    v3=np.cross(r11,r10)
    v3=v3/np.linalg.norm(v3)
    v4=np.cross(r10,r00)
    v4=v4/np.linalg.norm(v4)

    #pairwise dot products (cos of angle but already normalized)
    d1=np.dot(v1,v2)
    d1=1 if d1>1 else -1 if d1<-1 else d1
    d2=np.dot(v2,v3)
    d2=1 if d2>1 else -1 if d2<-1 else d2
    d3=np.dot(v3,v4)
    d3=1 if d3>1 else -1 if d3<-1 else d3
    d4=np.dot(v4,v1)
    d4=1 if d4>1 else -1 if d4<-1 else d4

    #pi/2 - (angle) 
    as1=np.arcsin(d1)
    as2=np.arcsin(d2)
    as3=np.arcsin(d3)
    as4=np.arcsin(d4)

    #Area(Qij) = arcsin of n1 * n2 + arcsin of n2 * n3, etc.
    alk=np.sign(np.dot(np.cross(ra,rb), r00))*(as1+as2+as3+as4)/(4*math.pi)
    return alk

#gauss_lk of a list rather than individual coordinates
def gauss_lkList(list):
    return gauss_lk(list[0], list[1], list[2], list[3])

#creates a random orthonormal 3*3 basis
def randomBasis(): 
    #generate random vector and normalize
    zv = np.random.normal(size=3)
    zv = zv/np.linalg.norm(zv) 
  
   #generate another random vector, find orthogonal projection, and normalize
    xv = np.random.normal(size=3)
    xv -= zv*np.dot(xv, zv) 
    xv = xv/np.linalg.norm(xv)
    
    #generate third unit vector orthogonal to first two
    yv = np.cross(zv, xv)

    print(np.array([xv, yv, zv]))
    return(np.array([xv,yv,zv]))

def vas_closed(verts):
    all_vas=[]
    
    #write as [[x,y,z],[x,y,z],...,[x,y,z]]
    v2=0
    nverts=len(verts)
    intersections=[]
    crossings=[]
    z=[]
    kcoords=[]
    
    r = randomBasis()
    #Takes vertices and puts into a new coordinate system based on a random vector. This gives us new vertices kcoords
    for k in range(0, nverts):
        xk=((r[0][0]*verts[k][0])+(r[0][1]*verts[k][1])+(r[0][2]*verts[k][2]))/(math.sqrt((r[0][0]*r[0][0])+(r[0][1]*r[0][1])+(r[0][2]*r[0][2])))
        yk=((r[1][0]*verts[k][0])+(r[1][1]*verts[k][1])+(r[1][2]*verts[k][2]))/(math.sqrt((r[1][0]*r[1][0])+(r[1][1]*r[1][1])+(r[1][2]*r[1][2])))
        zk=((r[2][0]*verts[k][0])+(r[2][1]*verts[k][1])+(r[2][2]*verts[k][2]))/(math.sqrt((r[2][0]*r[2][0])+(r[2][1]*r[2][1])+(r[2][2]*r[2][2])))
        kcoord=[xk,yk,zk]
        kcoords.append(kcoord)
    verts = kcoords           #these are the projected coordinates

    if nverts < 5:
        v2 = 0
        return v2
        
    for j in range(0, nverts): 
        vect01=np.array(verts[j])
        if j < nverts-1:
            vect02=np.array(verts[j+1])
        else:
            vect02=np.array(verts[0])
        
        for i in range(j+2, nverts):
            if i != j+nverts-1:
                vect11=np.array(verts[i])
                if i < nverts-1:
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

                #Cramer's rule to find intersection points in the projection
                A = np.array([ [x0, -x1], [y0, -y1] ])
                B = np.array([x11-x01, y11-y01])
                if np.linalg.det(A) != 0:
                    X = np.linalg.inv(A).dot(B)
                    if (0 <= X[0] <= 1) and (0 <= X[1] <= 1):       #Check if the potential intersection is within the edges
                        I = np.array([ X[0]*x0+x01, X[0]*y0+y01 ])
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

    verts.append(verts[0])
    for k in range(0,len(crossings)-1):
        for l in range(k+1,len(crossings)):
            u=crossings[k][0]
            v=crossings[k][1]
            w=crossings[l][0]
            x=crossings[l][1]

            #Check that crossings qualify based on 1 < 2, etc. and that they alternate correctly
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
                #Check for when one edge is the same for the two pairs
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

    verts.pop(len(verts)-1)
    # all_append(v2)
    
    # print(all_vas)
    v2 /= 2
    return v2

def vas_open(verts, niterations=100):
    nverts = len(verts)
    if nverts < 5:
        return 0        #If less than 5 points, the Vassiliev is 0

    sum_vas = 0
    # all_vas=[]

    for p in range(0, niterations):                 #n random projections of the knot so that the vas can be calculted n times and then averaged
        #write as [[x,y,z],[x,y,z],...,[x,y,z]]
        v2 = 0
        intersections, crossings, z, kcoords = [], [], [], []
        
        r = randomBasis()
        for k in range(0, nverts):
            xk = ((r[0][0]*verts[k][0])+(r[0][1]*verts[k][1])+(r[0][2]*verts[k][2])) / (math.sqrt((r[0][0]*r[0][0])+(r[0][1]*r[0][1])+(r[0][2]*r[0][2])))
            yk = ((r[1][0]*verts[k][0])+(r[1][1]*verts[k][1])+(r[1][2]*verts[k][2])) / (math.sqrt((r[1][0]*r[1][0])+(r[1][1]*r[1][1])+(r[1][2]*r[1][2])))
            zk = ((r[2][0]*verts[k][0])+(r[2][1]*verts[k][1])+(r[2][2]*verts[k][2])) / (math.sqrt((r[2][0]*r[2][0])+(r[2][1]*r[2][1])+(r[2][2]*r[2][2])))
            kcoords.append([xk,yk,zk])
        
        verts = kcoords           #these are the projected coordinates

        for j in range(0, nverts - 3):      #Happens nverts - 3 or around 1000 times
            vect01 = np.array(verts[j])
            vect02 = np.array(verts[j+1])

            for i in range(j+2, nverts-1):
                vect11 = np.array(verts[i])
                vect12 = np.array(verts[i+1])
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

                A=np.array([[x0,-x1],[y0,-y1]])     #Cramer's rule for intersection 
                if (np.linalg.det(A) != 0):
                    B = np.array([x11-x01,y11-y01])
                    X=np.linalg.inv(A).dot(B)

                    if (0 <= X[0] <= 1) and (0 <= X[1] <= 1):
                        I = np.array([X[0]*x0+x01, X[0]*y0+y01])
                        intersections.append(I)

                        zi0 = vect01[2]+(X[0]*(vect02[2]-vect01[2]))
                        zi1 = vect11[2]+(X[1]*(vect12[2]-vect11[2]))
                        
                        z.append([zi0,zi1])
                        crossings.append([j,i,X[0],X[1]]) 
        #print(f"\n{crossings}\n")\

        for k in range(0,len(crossings)-1):
            for l in range(k+1, len(crossings)):
                u=crossings[k][0]
                v=crossings[k][1]
                w=crossings[l][0]
                x=crossings[l][1]
                if (u<w) and (w<v) and (v<x) and ((z[k][0] > z[k][1] and z[l][0] < z[l][1]) or (z[k][0] < z[k][1] and z[l][0] > z[l][1])):
                    
                    v2 += float(np.sign(gauss_lk(np.array(verts[u]),np.array(verts[u+1]),np.array(verts[v]),np.array(verts[v+1])))) * float(np.sign(gauss_lk(np.array(verts[w]),np.array(verts[w+1]),np.array(verts[x]),np.array(verts[x+1]))))
                
                elif (z[k][0] > z[k][1] and z[l][0] < z[l][1]) or (z[k][0] < z[k][1] and z[l][0] > z[l][1]):
                        if (u==w) and (crossings[k][2]<crossings[l][2]) and (w<v) and (v<x):
                            #print("UW",crossings[k],crossings[l])
                            #print(u,v,w,x)
                            #print(verts[u])
                            #print(verts[u + 1])
                            #print(verts[v])
                            #print(verts[v + 1])
                            #print(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1]))))
                            v2 += float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]),np.array(verts[x + 1]))))
                        if (v==x) and (crossings[k][3]<crossings[l][3]) and (u<w) and (w<v):
                            #print("VX",crossings[k],crossings[l])
                            #print(u,v,w,x)
                            #print(verts[u])
                            #print(verts[u + 1])
                            #print(verts[v])
                            #print(verts[v + 1])
                            #print(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1]))))
                            v2 += float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]),np.array(verts[x + 1]))))
                        if (v==w) and (crossings[l][2]<crossings[k][3]) and (u<w) and (v<x):
                            #print("VW", crossings[k], crossings[l])
                            #print(u,v,w,x)
                            #print(verts[u])
                            #print(verts[u + 1])
                            #print(verts[v])
                            #print(verts[v + 1])
                            #print(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1]))))
                            #print(float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]), np.array(verts[x + 1])))))
                            v2 += float(np.sign(gauss_lk(np.array(verts[u]), np.array(verts[u + 1]), np.array(verts[v]), np.array(verts[v + 1])))) * float(np.sign(gauss_lk(np.array(verts[w]), np.array(verts[w + 1]), np.array(verts[x]), np.array(verts[x + 1]))))

        sum_vas += v2

        # all_append(v2)
        # print(v2)
    #end for p in range(0, niterations)

    #print(all_vas)
    vas_avg = (sum_vas/niterations) / 2
    return vas_avg

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
                r = randomBasis()
                pt.append(walk[i-1][0]+r[2][0])
                pt.append(walk[i-1][1]+r[2][1])
                pt.append(walk[i-1][2]+r[2][2])
                walk.append(pt)
            print(walk)
            vk = vas_open(walk)
            print(vk)
            vks.append(vk)
            vksum+=vk
        vkavg=vksum/len(vks)
        print(vks)
        print(vkavg)

#WII: Calculate second Vassiliev measure of knot and plot it
def plot_vas(knot, closed=False):
    vas = 0
    start_time = time.time()
    if closed:
        vas = vas_open(np.append(knot, knot[0]), 1)     #gives 0.5
    else:
        vas = vas_open(knot)
    end_time = time.time()

    print(f"The time elapsed was {end_time - start_time} seconds")
    print("\nThe Vassiliev Measure of the knot is " + str(vas) + "\n")

    x, y, z = [], [], []
    for i in range (0, len(knot)):
        x.append(knot[i][0])
        y.append(knot[i][1])
        z.append(knot[i][2])
    if closed:
        x.append(knot[0][0])
        y.append(knot[0][1])
        z.append(knot[0][2])

    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z)
    plt.show()

#WIII: calculate local second Vassiliev measure with varying number of atoms at once by calculating every length of the protein one index at a time
def vas_scan(protein):
    max_list = []
    interval = 100
    
    start_time = time.time()

    for scanlength in range(400, 601, 400):
        vas_list = []
        for start in range(0, len(protein) - scanlength - 1, interval):    #scan the protein at every possible starting index
            local_vas = vas_open(protein[start : start + scanlength])
            vas_list.append( local_vas )
            print(f"{local_vas} at {start}:{start + scanlength}")
        
        print(f"The Vas measures from 0 to {len(protein)} are {vas_list}")
        
        max_vas = max(np.abs(vas_list))                 #find max value of vas in given iteration, then record them
        max_start = vas_list.index(max_vas) * interval
        max_loc = [max_start, max_start + scanlength]
        print(f"The maximum Vassiliev for scanlength of {scanlength} is {max_vas}, at atoms {max_loc}\n")
        max_list.append([max_loc, max_vas])
    
    end_time = time.time()
    print(f"\nThe time elapsed was {end_time - start_time} seconds")
    return max_list

def plot_by_section(knot, section, interval):

    ax = plt.axes(projection='3d')
    for i in range(section[0], section[1] - interval, interval):
        x, y, z = [], [], []
        for j in range (i, i + interval, 1):
            x.append(knot[j][0])
            y.append(knot[j][1])
            z.append(knot[j][2])
        print(f"Knot from {i} to {i + interval} plotted")

        ax.plot3D(x, y, z)
        if (i + interval + interval > section[1]):
            for j in range (i + interval, section[1], 1):
                x.append(knot[j][0])
                y.append(knot[j][1])
                z.append(knot[j][2])
            print(f"Knot from {i + interval} to {section[1]} plotted")
    
    plt.show()

##########################

# trefoil stuff
trefoil = [[1, 0, 0],
            [4, 0, 0],
            [1, 6, 2],
            [0, 2, -5],
            [5, 2, 5],
            [4, 6, -2]]

plot_vas(trefoil, True)
#vas is 1.0 for closed trefoil

    # trefoil = [[1, 0, 0],
    #             [4, 0, 0],
    #             [1, 6, 2],
    #             [0, 2, -5],
    #             [5, 2, 5],
    #             [4, 6, -2],
    #             [0.5, 0.5, 0.5]]

    # plot_vas(trefoil)
    # #vas is 0.995 for open trefoil

# proteinName = ["6zge"]
# proteinDF = pd.read_csv(fr'Coordinates\{proteinName}.csv')
# proteinList = proteinDF.values.tolist()                 #change df to a list of atoms' coordinates
# print(f"Protein has {len(proteinList)} CA atoms\n")

# start_time = time.time()            #start
# spike_vas = vas_open(proteinList)
# end_time = time.time()              #end
# print(f"\nThe time elapsed for the vas_open method was {end_time - start_time} seconds")
# pd.Dataframe([proteinName, len(proteinList), (end_time - start_time)])  #try to create a DF out of all these things

# vas_list = vas_scan(spikeList)
# print(f"vas list: {vas_list}")

# plot_by_section(spikeList, [0, len(spikeList)], 200)
# interval = 200
# for i in range(0, len(spikeList), interval):
#     plot_by_section(spikeList, [i, i + interval], interval)         #change 0 to i to have it show individually
