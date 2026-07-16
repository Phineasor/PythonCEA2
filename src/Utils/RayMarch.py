"""
@author: phineas
"""
#Important imports
#External Libraries
import math as m
import numpy as np
import matplotlib.pyplot as plt

#Internal libraries
from .EngineGeometry import RatL, LT
from ..Inputs import InputValues as IV
from .Bisect import Bisect



#This fuiction gets the ray data from an arbitrary location across the engine in terms of cell number, and an arbitrary pair of angles measures from pointing at the top of the engine
#the other directionality is irrelivent because the engine is symetric from any plane around the axis
#We will obtain the power values at linspace location across the engine section, however we will quite probably need to polynomialy fit data to get some points inbetween

def getRay(x, theta1, theta2, test = False):
    #determines the distance and the radial distance from the z axis at each locatio. Y will be straight twards the wall in a positive direction, X is measured from the injector face in a posive manor following into the end
    AxialDistances = np.linspace(0, LT, IV.CellNum)
    RadiusVal = [0.0]*IV.CellNum
    for i in range(IV.CellNum):
        RadiusVal[i] = RatL(AxialDistances[i])

    #computes the thickness of a "mesh cell"
    Lx = LT/(IV.CellNum-1)

    #Next we need to obtain the vector normal to the plane that we are working with for the angle of the chamber wall. This is found with the cross product

    dx = 10**(-6) #value for finding the 2d derivative in the xy plane
    
    #will use the slope of the point farther down the engine unless is the last cell, then it will be the one behiend it
    if not(x == (IV.CellNum-1)):
        p1 = np.array([AxialDistances[x], RadiusVal[x], 0])
        p2 = np.array([AxialDistances[x]+dx, RatL(AxialDistances[x]+dx), 0])
        
        Slope = p2-p1
        NormalVector = np.cross(Slope, np.array([0, 0, 1]))
        UNormalV = NormalVector/np.linalg.norm(NormalVector)
    else: 
        p1 = np.array([AxialDistances[x], RadiusVal[x], 0])
        p2 = np.array([AxialDistances[x]-dx, RatL(AxialDistances[x]-dx), 0])
    
        Slope = p1-p2
        NormalVector = np.cross(Slope, np.array([0, 0, 1]))
        UNormalV = NormalVector/np.linalg.norm(NormalVector)
    #p1 and p2 are flipped to produce the negetive NormalVector such that the vectors used to construct it can be used in a change of basis matrix
    #we now have the unit normal vector that is orthogonal to the wall at the locaiton, this is used to cefine the zy plane reziding in this region
    
    #we now need the 3d change of basis matrix so that we can transform these vectors vetween the engine centered frame, and the local wall frame, this is important for doing the rotation matrix on the vector to point in a new direction.
    #conviniently we know all three basis vectors, the normal vector for the x axis, adn the two vectors used to construct it in the array
    USlope = Slope/np.linalg.norm(Slope)
    W = np.array([
        [(UNormalV[0]), (USlope[0]), (0)],
        [(UNormalV[1]), (USlope[1]), (0)],
        [(UNormalV[2]), (USlope[2]), (1)]
    ])
    #now in theory this is the transformation matrix so lets normalize the vector that we need to rotate
    RotateVector = np.linalg.inv(W) @ UNormalV
    
    
    #now we need to apply two rotation matricies to rotate this vector dependent on the theta1 and theta2
    RmatPitch   = np.array([
        [(m.cos(-theta1)), (-m.sin(-theta1)), (0)],
        [(m.sin(-theta1)), (m.cos(-theta1)), (0)],
        [(0), (0), (1)]
    ])
    RmatRoll = np.array([
        [(1), (0), (0)],
        [(0), (m.cos(theta2)), (-m.sin(theta2))],
        [(0), (m.sin(theta2)), (m.cos(theta2))]
    ])
    #the theta1s in RmatYaw are to ensure that for 0-90 theta2 the value remains positive, this is not a requirement but it makes it more convinient
    
    #now we just compute the fully rotated vector
    RotatedVector = RotateVector @ RmatPitch @ RmatRoll #get rotated idiot
    RotatedVectorUnrot = W @ RotatedVector 
    
    line = lambda t, point, slope:[(point[0]+t*slope[0]), (point[1]+t*slope[1]), (point[2]+t*slope[2])]

    tval1 = Bisect(getPoint, (10**(-6)), 10**10, (10**(-5)), p1, RotatedVectorUnrot, False)
    tval2 = Bisect(getPoint, (10**(-6)), 10**10, (10**(-5)), p1, RotatedVectorUnrot, True)
    #print(tval1)
    #print(tval2)
    if (tval1 < tval2):
        tval = tval1
    else:
        tval = tval2

    intersect = np.array(line(tval, p1, RotatedVectorUnrot))
   
    #handles the case where the ray leaves the confines of the engine.
    if intersect[0] < 0:
        scale = ((-p1[0])/RotatedVectorUnrot[0])
        intersect = line(scale, p1, RotatedVectorUnrot)
    if intersect[0] > LT:
        scale = ((LT-p1[0])/RotatedVectorUnrot[0])
        intersect = line(scale, p1, RotatedVectorUnrot)
    
 
    
    #Determins all of the array points that the ray passes through, returns position and point in the array
    locations = np.array([p1])
    arraynums = np.array([x])
    point = p1[0]+Lx/4 #adds some small length based amount to ensure no equivilence garbage happens
    
    i=0
    while not (((Lx*(i+1)+point)>intersect[0])):
        #print(i)
        #Will scan the array of distances to find what number the current section of this ray has passed through
        j = 0
        Found = False
        while (not Found):
            if ((AxialDistances[j] > (Lx*(i)+point)) and (AxialDistances[j] < (Lx*(i+1)+point))):
                Found = True  
            else:
                j+=1
            #print('p1: ' + str((Lx*(i)+point)) + ' | p2: ' + str((Lx*(i+1)+point)))
        #now j is the value in the axial distance array that the ray is next to pass through, so we add it to the values that we pass through
        #print('j: ' + str(j))
        arraynums = np.append(arraynums, np.array([j]))
        
        #we now wish to use this infromation to determine the precise location of this data point as well in 3D space
        scale = abs(p1[0]-AxialDistances[j])/RotatedVectorUnrot[0]
        locations = np.append(locations, np.array([line(scale, p1, RotatedVectorUnrot)]), axis = 0)
        i+=1
    #Finds the last element, this ensures that it is always included, even for straigt rays.
    j = 0
    Found = False
    while (not Found):
            if ((AxialDistances[j] > (intersect[0]-Lx*0.6)) and (AxialDistances[j] < (intersect[0]+Lx*0.6))):
                Found = True  
            else:
                j+=1
    arraynums = np.append(arraynums, np.array([j]))
    locations = np.append(locations, np.array([intersect]), axis = 0)
    
    
    
    ray = [x, theta1, theta2, np.flip(locations, axis = 0), np.flip(arraynums, axis = 0)]
    if test:
        return [p1, RotatedVectorUnrot, intersect, (intersect[1]**2+intersect[2]**2)**0.5, tval, i, abs(i-x)]
    else:
        return ray



#Function that determines the difference between a line at some paramaterized point t, and some 
def getPoint(t, point, slope, negQM):
    if negQM:
        point = [point[0], -point[1], point[2]]
        slope = [slope[0], -slope[1], slope[2]]
    pointBeingChecked = [(point[0]+t*slope[0]), (point[1]+t*slope[1]), (point[2]+t*slope[2])]

    RootingTerm = ((RatL(pointBeingChecked[0]))**2-(pointBeingChecked[2])**2)
    #Prevents a negetive root, and provides the needed info to let the bisect function know that it has overshot to keep going.
    if((RootingTerm < 0) or (RootingTerm > RatL(0)**2+1)):
        return 10

    #sometimes the line will intersect on the same side of the engine so we need a positive and negetive version fo this math.
    ChamberZ =  m.sqrt(RootingTerm)

    #print("RT: " + str(RootingTerm))
    #print("point: " + str(pointBeingChecked[1]))
    #print("chamber: " + str(ChamberZ))
    #print("val: " + str((pointBeingChecked[1]-ChamberZ)))
    return (pointBeingChecked[1]-ChamberZ)


#print(getRay(0, (10*(m.pi/180)), (20*(m.pi/180))))
#print(getRay(0, (89*(m.pi/180)), (0*(m.pi/180))))
#print(getRay(0, (85*(m.pi/180)), (0*(m.pi/180))))
#print(getRay(0, (-85*(m.pi/180)), (0*(m.pi/180))))
#print(getRay(0, (-89.9999*(m.pi/180)), (0*(m.pi/180))))
#print(getRay(0, (85*(m.pi/180)), (20*(m.pi/180))))
#print(getRay(245, (85*(m.pi/180)), (85*(m.pi/180))))
#print(getRay(0, (75*(m.pi/180)), (0*(m.pi/180)))[-2])
#print(getRay(0, 0, 0))


test = False
#print(LT)

if test:
    #RAYMARCH testing 
    plt.rcParams['figure.dpi'] = 500
    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111,projection='3d')
    ax.set_facecolor('black')


    num = 100
    Edist = np.linspace(0, LT, num)
    i = 0
    while(i < 360):
        j = 0
        list1 = [0]*num
        list2 = [0]*num
        list3 = [0]*num
        while(j < num):
            list1[j] = Edist[j]
            list2[j] = m.sin((m.pi/180)*i)*RatL(Edist[j])
            list3[j] = m.cos((m.pi/180)*i)*RatL(Edist[j])
            j += 1
        ax.plot(list1,list2,list3, color='white',linestyle='-',linewidth=1) 
        i += 8

    
    #testray = getRay(0, (10*(m.pi/180)), (20*(m.pi/180)))
    #testray = getRay(0, (85*(m.pi/180)), (20*(m.pi/180)))
    #testray = getRay(0, (-89.99999*(m.pi/180)), (0*(m.pi/180)))
    #testray = getRay(0, (0*(m.pi/180)), (0*(m.pi/180)), True)
    #testray = getRay(245, (85*(m.pi/180)), (85*(m.pi/180)), True)
    #testray = getRay(0, (89*(m.pi/180)), (0*(m.pi/180)), True)
    #print(testray)

    ax.scatter(testray[0][0], testray[0][1], testray[0][2], color='green',linestyle='--', linewidth=0.2)
    ax.scatter(testray[2][0], testray[2][1], testray[2][2], color='red',linestyle='--', linewidth=0.2)

    line = lambda t, point, slope:[(point[0]+t*slope[0]), (point[1]+t*slope[1]), (point[2]+t*slope[2])]


    nums = np.linspace(0, 10, 100)
    linex = [0]*100
    liney = [0]*100
    linez = [0]*100

    for i, val in enumerate(linex):
        num = line(nums[i], testray[0], testray[1])
        linex[i] = num[0]
        liney[i] = num[1]
        linez[i] = num[2]

    ax.plot(linex, liney, linez, color='blue',linestyle='-', linewidth=0.5)

    test2 = False
    
    if test2:
        testray2 = getRay(0, (75*(m.pi/180)), (0*(m.pi/180)), False)
        for i in testray2[3]:
            ax.scatter(i[0], i[1], i[2], color='yellow', marker ='.', s = 0.1)
    


    ax.set_xlim3d(-LT/2, LT/2)
    ax.set_ylim3d(-LT/2, LT/2)
    ax.set_zlim3d(-LT/2, LT/2)
    ax.grid(False)
    plt.title("3D engine model")
    plt.legend(fancybox=False, shadow=True, framealpha=1,fontsize='small',loc='lower left')
    plt.show()



'''
count = 0
ihatethis = 0
j = -89
k = 0
while ihatethis < 250:
    j = -89
    while j < 90:
        k = 0
        while k < 90:
            test = getRay(ihatethis, (j*(m.pi/180)), (k*(m.pi/180)))
            k+=2
            count += 1
            if (count % 1000 == 0):
                print(count)
        j+=2
    ihatethis+=1    
'''