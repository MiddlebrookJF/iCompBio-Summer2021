#THIS IS MEANT FOR USE IN THE CLUSTER

import os
import math
import time
import numpy as np
import pandas as pd
import glob
from multiprocessing import Pool
from functools import partial

#finds scalar triple product of 3 3-D vectors
def tripleProduct(a,b,c): 
    return np.dot(np.cross(a,b),c)

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

    return [xv,yv,zv]

#calculates crossings given walk and two points (vectors to next points from those points)
def crossing (walk, i, k):
    #set up variables, separating out coordinates
    x00, x01, x10, x11 = [walk[i][0], walk[i+1][0], walk[k][0], walk[k+1][0]]
    y00, y01, y10, y11 = [walk[i][1], walk[i+1][1], walk[k][1], walk[k+1][1]]
    z00, z01, z10, z11 = [walk[i][2], walk[i+1][2], walk[k][2], walk[k+1][2]]
    
    #Cramer's rule to find parametrized intersection
    A = [[x01-x00, x10-x11],[y01-y00, y10-y11]]
    B = [[x10-x00],[y10-y00]]
    if np.linalg.det(A) != 0:
        paramInt = np.matmul(np.linalg.inv(A),B) 
        frac1 = paramInt[0] 
        frac2 = paramInt[1]
        if (0 <= frac1 <= 1) and (0 <= frac2 <= 1): #checks if intersection is within vectors
            z0 = z00 + frac1*(z01-z00)
            z1 = z10 + frac2*(z11-z10)
            return ([i, k, i if z0>z1 else k, frac1, frac2]) #return two original points, the overstrand, and the fractions
        else:
            return None
    else:
        return None

#checks each of the possible conditions for a pair of crossings, returns true if any are true
def vas_conditions(cross1, cross2):
    conditions = []
    i,j,k,l = [cross1[0],cross2[0],cross1[1],cross2[1]]

    conditions.append(i<j<k<l)
    conditions.append(i==j<k<l and cross1[3]<cross2[3])
    conditions.append(i<j==k<l and cross2[3]<cross1[4])
    conditions.append(i<j<k==l and cross1[4]<cross2[4])
    
    return any(conditions)
    

#calculates second Vassiliev measure of just one projection, using crossing and vas_conditions
def vas_proj(walk, proj=[[1,0,0],[0,1,0],[0,0,1]]):
    nverts = len(walk)
    vas_sum = 0
    IList = []
    
    #transform the coordinates of the walk
    walk = np.matmul(walk, np.transpose(proj).tolist())
    
    #create list of pairs to check
    pairs=[]
    for j in range(nverts-3):
        for k in range(j+2,nverts-1):
            intersection = crossing(walk,j,k)
            if(intersection != None):
                IList.append(intersection)

    #for each pair of crossings, check if ordered correctly 
    for cross1 in IList: 
        for cross2 in IList: 
            i,j,k,l = [cross1[0],cross2[0],cross1[1],cross2[1]] #i,j,k,l are indices of walk with a crossing
            if (vas_conditions(cross1, cross2)): #ensure coordinates are in order
                if (cross1[2]==i and cross2[2]==l) or (cross1[2]==k and cross2[2]==j): #one is under and one is over
                    signedprod1 = np.sign(tripleProduct(walk[i+1]-walk[i], walk[k+1]-walk[k],walk[i]-walk[k]))
                    signedprod2 = np.sign(tripleProduct(walk[j+1]-walk[j], walk[l+1]-walk[l],walk[j]-walk[l]))
                    vas_sum += signedprod1*signedprod2 #add the signed product (+/- 1)
    
    return(vas_sum * 0.5)

#estimates second Vassiliev measure of open chain, using vas_proj, crossings, vas_conditions
#takes a walk as parameter, number of trials and list of random bases are optional
def vas_open(chain, trials=1000, size=50, poolNum=4):
    random_list = []
    vas_list = []
    if len(chain) < 4: return 0
    for i in range(trials):
        random_list.append(randomBasis())

    for j in random_list:
        vas_list.append(vas_proj(chain, j))
    
    #print(vas_list)
    vas_sum = sum(vas_list)/trials
    return vas_sum

def vas_open_parallel(chain, trials=1000, size=50, poolNum=4):
    random_list = []
    if len(chain) < 4: return 0
    
    for i in range(trials):
        random_list.append(randomBasis())

    if __name__== '__main__':
        #iterate over the list with multiple processors
        p = Pool(poolNum)
        part = partial(vas_proj, chain)
        result = p.map(func = part, iterable = random_list, chunksize = size)
        if(result != None):
            #print(result)
            vas_sum = sum(result)/trials
        p.close()
        p.join()
        return vas_sum

#vas of either open or closed chain, passed as boolean; default is open
def vas_measure(walk, closed = False, trials=1000):
    if (closed):
        walk.append(walk[0])
        return vas_open_parallel(walk, trials=1)
    return vas_open_parallel(walk, trials=trials)

#calculate runtime
def runtime (startTime):
   return time.time()-startTime

#WIII: calculate local second Vassiliev measure with varying number of atoms at once by calculating every length of the protein one index at a time
def vas_scan(protein, numProjections=1000):
    max_list = []
    interval = 50
    pLength = len(protein)

    for scanlength in range(200, 601, 200):
        sTime = time.time()
        vas_list = []

        for start in range(0, pLength - scanlength, interval):     #scan the protein in range of length scanlength starting at 'start' which += by interval
            upperbound = start + scanlength
            local_vas = vas_open_parallel(protein[start : upperbound], trials=numProjections)
            vas_list.append( local_vas )
            print("{local_vas} at {start}:{upperbound}".format(local_vas=local_vas, start=start, upperbound=upperbound))

            if(upperbound + interval > pLength):                   #last iteration
                if(pLength < upperbound + interval - (interval/2)):           #Determines range of last scan
                    start = pLength - interval
                else: start += interval
                
                local_vas = vas_open(protein[start : pLength])
                vas_list.append( local_vas )
                print("{local_vas} at {start}:{pLength}".format(local_vas=local_vas, start=start, pLength=pLength))

        print("The Vas measures from 0 to {pLength} are {vas_list}".format(pLength=pLength, vas_list=vas_list))
        
        max_vas = max(np.abs(vas_list))                 #find max value of vas in given iteration, then record them
        max_start = vas_list.index(max(vas_list)) * interval
        max_loc = [max_start, max_start + scanlength]
        print("The maximum Vas for scanlength {scanlength} is {max_vas}, at atoms {max_loc}".format(scanlength=scanlength, max_vas=max_vas, max_loc=max_loc))
        max_list.append([max_loc, max_vas])

        print("{execMin} minutes of runtime for {scanlength} scanlength\n".format(execMin=runtime(sTime)/60, scanlength=scanlength))
    
    return max_list

##########################

### Proteins Vas_Scan FOR CLUSTER ###

proteins = glob.glob('Coordinates/*.csv')
numProjections=1000
for proteinPath in proteins[int(len(proteins)/3) : int(len(proteins)*2/3)]:
    proteinDF = pd.read_csv(proteinPath)
    proteinName = proteinPath[-8:-4]
    proteinList = proteinDF.values.tolist()                 #change df to a list of atoms' coordinates

    ## Overall Vas ##
    startTime = time.time()
    value = vas_open_parallel(proteinList, trials=numProjections)
    execTime = runtime(startTime)
    if(value != None):
        thisDF=pd.DataFrame([[proteinName, len(proteinList), value, execTime, numProjections]],
            columns=['Name','NumCaAtoms','Vassiliev','RuntimeSeconds','NumProjections'])
        with open(f"Vas-Data/NewVas{numProjections}.csv", mode='a') as f:
            thisDF.to_csv(f, header=f.tell()==0, index=False)
