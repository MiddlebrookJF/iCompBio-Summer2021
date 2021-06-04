#import statements
import os
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d
import Read_File as read
import time
from multiprocessing import Pool, Process
from functools import partial



#FUNCTION DEFINITIONS

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

   #calculate and return final value
    alk = np.sign(tripleProduct(ra,rb,r00))*(as1+as2+as3+as4)/(4*math.pi)
    #sgn=np.sign(aux)
    return alk

#gauss_lk of a list rather than individual coordinates
def gauss_lkList(list):
    return gauss_lk(list[0], list[1], list[2], list[3])

#finds writhe of an entire chain; input is a chain of coordinates as an np.array
#uses gauss_lk and gauss_lk_list
def compute_wr(walk, size = 1000, poolNum = 2):
    wr=0.0
    acn=0.0
    nvertices=len(walk)
    
    #multiprocessing unit
    # if __name__ == '__main__':
        #create list with every set of vertices to pass through gauss_lk
    pairlist=[]
    for j in range (nvertices-3):
            u1=walk[j]
            u2=walk[j+1]
            for i in range (j+2,nvertices-1):
                v1=walk[i]
                v2=walk[i+1]
                pairlist.append([u1,u2,v1,v2])
    
    #iterate over the list with multiple processors
    p = Pool(poolNum)
    result = p.map(func = gauss_lkList, iterable = pairlist, chunksize = size)
    if(result != None):
        wr = 2*sum(result)
        acn = 2*sum(abs(i) for i in result)
    p.close()
    p.join()
    return wr, acn

#creates a random orthonormal 3*3 basis
def randomBasis(): 
    #generate random vector and normalize
    zv = np.random.normal(size=3)
    zv /= np.linalg.norm(zv) 
  
   #generate another random vector, find orthogonal projection, and normalize
    xv = np.random.normal(size=3)
    xv -= zv*np.dot(xv, zv) 
    xv /= np.linalg.norm(xv)
    
    #generate third unit vector orthogonal to first two
    yv = np.cross(zv, xv)

    return(np.array([xv,yv,zv]))
           
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
def vas_proj(walk, proj=np.array([[1,0,0],[0,1,0],[0,0,1]])):
    nverts = len(walk)
    vas_sum = 0
    IList = []
    
    #transform the coordinates of the walk
    walk = np.matmul(walk, np.transpose(proj))
    
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
def vas_open(chain, trials=100, size=10, poolNum=2):
    random_list = []
    vas_list = []
    for i in range(trials):
        random_list.append(randomBasis())

    for j in random_list:
        vas_list.append(vas_proj(chain, j))
    
    #print(vas_list)
    vas_sum = sum(vas_list)/trials
    return vas_sum

def vas_open_parallel(chain, trials=100, size=10, poolNum=2):
    random_list = []
        
    for i in range(trials):
        random_list.append(randomBasis())

    # if __name__== '__main__':
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
def vas_measure(walk, closed = False):
    if (closed):
        return vas_open(np.append(walk, walk[0]), trials=1)
    return vas_open(walk)

#creates a random walk of length n, each segment length 1
def randwalk(n):
    walk=[[0,0,0]]
        #fig = plt.figure()
        #ax = plt.axes(projection='3d')
    for i in range (n):
        r=randomBasis()
        walk.append(walk[i]+r[2])
    #print(walk)
    walk = np.array(walk)
    return walk

#vas of number of random walks of length n
def vas_randwalk(n, sample_size=10):
    vk_list = []
    vksum = 0
    
    #generates sample_size random walks
    for t in range(sample_size):
        walk_sample = randwalk(n)
        vk = vas_open(walk_sample)
        vk_list.append(vk)
        vksum+=vk
    vkavg=vksum/sample_size
    print('Vassiliev measure of each walk: \n', vk_list)
    print('Average Vassiliev measure:', vkavg)

#plot chain
def plot3D(walk):
    x=[]
    y=[]
    z=[] #initializes x,y, and z; leave empty
    for i in range (len(walk)):
        x.append(walk[i][0])
        y.append(walk[i][1])
        z.append(walk[i][2])
    fig=plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z)
    plt.show()

#calculate runtime
def runtime (startTime):
   return time.time()-startTime

#### RUN COMPUTATIONS BELOW HERE

# proteins = read.readAll_pdb(fr"Proteins-PDB") #insert file path of pdb files here (/*folder/*folder/etc), or leave blank
# proteins = sorted(proteins.items(), key=lambda x: len(x[1])) #returns tuples of (protein, chain) sorted by chain length

# #estimate vas of each point
# vasValues = {}
# for tuple in proteins:
#     # value = None
#     # while value == None: #ensure an actual return result before moving on
#     startTime = time.time()
#     value = vas_open_parallel(tuple[1],1000, size=20)
#     execTime = runtime(startTime)
#     if(value!=None):
#         print (tuple[0], ':' , len(tuple[1]))
#         print('Vas: %f' %(value))
#         print('Runtime: %f \n'%(execTime))
#         vasValues.update({tuple[0]:value})

##match protein name with point, append to plotted tuple
#foldData=["1qpu",5.30,"5mbn",4.83,"1lmb",4.78,"2pdd",4.20,"1hrc",3.80,"1imq",3.16,"2abd",2.85,"2vil",3.25,"2hbb",2.87,"1ubq",3.19,"1cis",1.75,"1urn",2.53,"3gb1",2.46,"2ptl",1.78,"1fkf",0.60,"1hdn",1.17,"1afi",0.26,"1aps",-0.64,"1csp",2.84,"1tit",1.51,"1shf",1.97]
#plotpoints = []
#foldingRates = {}
#for i in range(int(len(foldData)/2)):
#    foldingRates.update({foldData[2*i]:foldData[2*i+1]})
#for pair in proteins:
#    plotpoints.append((pair[0], vasValues[pair[0]]), foldingRates[pair[0]] )

##create scatterplot
#x=[]
#y=[]#initializes x,y

#for i in plotpoints:
#    x.append(i[1])
#    y.append(i[2])
#    plt.annotate(i[0], (i[1],i[2]))
#plt.scatter(x, y)
#plt.title('Comparing Vas and Folding rates of simple proteins')
#plt.xlabel('Vas of chain')
#plt.ylabel('Folding rate')
#plt.show()