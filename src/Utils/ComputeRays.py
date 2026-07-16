# File to generate radiative power values
# fmt: off
"""
@author: phineas
"""

#External modules
import numpy as np
import math as m
from multiprocessing import shared_memory

#Internal modules
from cea import *
import InputValues as IV
from RayMarch import getRay
from concurrent.futures import ProcessPoolExecutor

#some universal constants
c = 299792458 #speed of light in meters per second
h = 6.62607015*10**(-34) #planks constant
kb = 1.380649*10**(-23)

ceaOut = runCEA()
val = AxialValues(ceaOut[0].T, ceaOut[0].P, ceaOut[0].density, ceaOut[0])

Test = False
#Alters memory of retrieving the AbsCoefAray for memory maping
if IV.Memmap and not Test:
    AbsCoefAray = np.load(IV.AbsCoefName, allow_pickle=True, mmap_mode = 'r')
elif not IV.Memmap and not Test:
    existing_shm = shared_memory.SharedMemory(name = 'AbsCoefDataMemoryBuffer')
    AbsCoefAray = np.ndarray((IV.CellNum, ), dtype = np.float64, buffer = existing_shm.buf)
else:
    AbsCoefAray = np.load(IV.AbsCoefName, allow_pickle=True)
length = len(AbsCoefAray[0][0])

#returnes the integrated value for one ray t
def CompRay(x, theta1, theta2):
    x = int(x)
    #print("test: " + str([x, theta1, theta2]))
    if (x == 0) and (theta1 < 0):
        return [x, theta1, theta2, 0]
    elif (x == IV.CellNum) and (theta1 > 0):
        return [x, theta1, theta2, 0]
    else:
        RadiativePower = 0
        Ray = getRay(x, (theta1*(m.pi/180)), (theta2*(m.pi/180)))
        #runs the needed number of times to produes the 
        i = 0
        p1 = CompRayAtWavenumber(i, [Ray[3], Ray[4]])
        wavelength1 = (1/AbsCoefAray[0][0][i])/100
        while i < length-1: #Integrate over wavenumber
            p2 = CompRayAtWavenumber(i+1, [Ray[3], Ray[4]])
            wavelength2 = (1/AbsCoefAray[0][0][i+1])/100
        
            RadiativePower += (wavelength1-wavelength2)*(p1+p2)/2
            
            p1 = p2
            wavelength1 = wavelength2
            if i%500000 == 0:
                print("compray i: " + str(i))
            i+=1
            
    print("Ray_Done----------------------------------------------------")
    return [x, theta1, theta2, RadiativePower] #[x, theta1, theta2, watts]

def CompRayAtWavenumber(wavenumberIndex, Ray):
    Intensity = 0
    wavenumber = AbsCoefAray[0][0][wavenumberIndex]
    i = 0
    p1 = plank(val[0][i], (1/wavenumber)/100)*100*AbsCoefAray[Ray[1][i]][1][wavenumberIndex]*m.exp(-opticalDpeth(wavenumberIndex, Ray, i))
    while i < len(Ray[0])-1:
        p2 = plank(val[0][i+1], (1/wavenumber)/100)*100*AbsCoefAray[Ray[1][i+1]][1][wavenumberIndex]*m.exp(-opticalDpeth(wavenumberIndex, Ray, i+1))     
        
        distance = (2.54*((Ray[0][i][0]-Ray[0][i+1][0])**2+(Ray[0][i][1]-Ray[0][i+1][1])**2+(Ray[0][i][2]-Ray[0][i+1][2])**2)**0.5)/100
        Intensity += ((p1+p2)/2)*distance
        i+=1
        p1 = p2
    return Intensity



#planks function dependent on temp in kelvin, and wavelength in meters
def plank(T, y):
    t1 = (2*h*c**2)/(y**5)
    exp1 = h*c
    exp2 = y*kb*T
    
    t2 = 1/(math.exp(exp1/exp2)-1)
    return t1*t2



#Function for computing the optican depth form 1 array points on the ray values it always assumes its to the end and some wavenumber index\
#this functionaly integrates absorbances from s to the end of the array, be it 
def opticalDpeth(wavenumberIndex, Ray, s): #ray should be [locations, arraynums]
    OD = 0
    i = s
    p1 = AbsCoefAray[Ray[1][i]][1][wavenumberIndex]
    while i < len(Ray[0])-1:
        p2 = AbsCoefAray[Ray[1][i+1]][1][wavenumberIndex]
            
        distance = 2.54*((Ray[0][i][0]-Ray[0][i+1][0])**2+(Ray[0][i][1]-Ray[0][i+1][1])**2+(Ray[0][i][2]-Ray[0][i+1][2])**2)**0.5 #conversion for in to cm at the beginning, absorbance is in cm
        OD += distance*(p1+p2)/2
        i+=1
        p1 = p2
    return OD


#testray = getRay(0, (75*(m.pi/180)), (0*(m.pi/180)))
#print(AbsCoefAray[testray[4][0]][0][200000])
#print(AbsCoefAray[testray[4][0]][1][200000])
#print(opticalDpeth(200000, [testray[3], testray[4]], 5)) 

#print(CompRayAtWavenumber(800000, [testray[3], testray[4]]))
#print(CompRay(0, (75*(m.pi/180)), (0*(m.pi/180))))
#print(CompRay(0, (75*(m.pi/180)), (0*(m.pi/180))))
