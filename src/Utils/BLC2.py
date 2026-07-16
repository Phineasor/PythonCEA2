#Good luck, I sure wish I had some
#This code takes the form of a function that can be called at any time, it determines the temperature profile of the wall along with whatever else it outputs (I dont know what that may be while writing this)
"""
@author: phineas
"""
#Important imports
#External Libraries
import math as m
import numpy

#Internal libraries
from ..Utils import Injector as Inj
from ..Inputs import InputValues as IV
from .EngineGeometry import RatL, Dt, LT, Lt
from .Bisect import Bisect
from .RayMarch import getRay
from .cea import *

#Random constants
in2m = 0.0254
d2r = m.pi/180


#This function takes in a location along the engine, and engine conditions and returns a radiative heat flux value using some... interesting methodologies.
def getEmittance(xi, ray):
    #We will need a to determing the emissivity from every ray from any angle in the engine with a seemingly unique aproach Grissom mentions an ideal methodology, we will use something simmilar
    #Averaging values of emissivity from various angles, we will just use the mainstream combustion gas temperature for now, any values in here that would have less than this value of temperature
    #would need to be emitted through a large volume of gas causing much of it to get absorbed or scattered.
    
    #we need to loop through all optical path lengths to    
    ray = getray(xi, theta1, theta2)
     
    
    return ray

#we must determine the "optical path length" in bar*cm, first we get some ray info at an angle this function determines the emittance 
#froma  single ray direction
def getEmittanceAtRay(x1, theta1, theta2, gas):
    #needed values of the gas propties thoughout the engine
    val = AxialValues(gas.T, gas.P, gas.density, gas)
    
    #we must determine the "optical path length" in bar*cm, first we get some ray info at an angle
    ray = getRay(x1, theta1, theta2)
    
    #we need to integrate partial pressure in bar along the length of this array to get the optical path length
    #There will be two pLs at the end for both gasses given they are partial pressures
    p1w = val[1][ray[4][0]]*gas.X[gas.species_index('H2O')] #starting pressure value at the beginning of the trapizoidal sum
    p1c = val[1][ray[4][0]]*gas.X[gas.species_index('CO2')]
    pLw = 0 #starting value for optical path length
    pLc = 0
    #trapizoidal integration to find optical path length in bar*cm
    for i in range(len(ray[3])-1):
        p2w = val[1][ray[4][i+1]]*gas.X[gas.species_index('H2O')] #pressure value at next point along ray
        p2c = val[1][ray[4][i+1]]*gas.X[gas.species_index('CO2')]
        pLpw = 2.54*(1/100000)*cartDist(ray[3][i], ray[3][i+1])*(p1w+p2w)*0.5 #constants to convert to cm and bar from in and pa
        pLpc = 2.54*(1/100000)*cartDist(ray[3][i], ray[3][i+1])*(p1c+p2c)*0.5
        
        #add to the term, and switch to the next term
        pLw += pLpw
        pLc += pLpc
        p1w = p2w
        p1c = p2c

    #PL should not be an optical path length in bar cm
    #now we must use this optical path length to determine the emissivity according to leckner using several polynomial fits
    #we assume that most of the emission happens in the chamber, so we use the chamber temperature for these values of optical path length
    #This is probably correct from the simmilar values between the optical path length in only the chamber, and straight through the whole nozzle, especially for near square sized combustion chambers like ours
    #future comparisons will be made with an average temperature along the ray aswell.
    T = gas.T
    
    #------------------------------------------------------------
    #T = 1001
    #pLw = 200
    #------------------------------------------------------------
    
    #For these correlations we must use two substetuded variables tao and lambda
    tao = T/1000
    lamw = m.log10(pLw) #not too sure if this is ln or not
    lamc = m.log10(pLc)
    
    #not we must apply the polynomial equations the coefficients are given, these coefficients are all listed as arrays here
    waterCs1 = [[-2.2118, -1.1987, 0.035596],
                [0.085667, 0.93048, -0.14391],
                [-0.10838, -0.17156, 0.045915]]
    
    CO2Cs1 = [[-3.9893, 2.7669, -2.1081, 0.39163],
              [1.2701, -1.1090, 1.0195, -0.21897],
              [-0.23678, 0.19731, -0.19544, 0.044644]]
    
    CO2Cs2 = [[-3.7981, 2.7353, -1.9882, 0.31054, 0.015719],
              [1.9326, -3.5932, 3.7247, -1.4535, 0.20132],
              [-0.35366, 0.1766, -0.84207, 0.39859, -0.063356],
              [-0.080181, 0.31466, -0.19973, 0.046532, -0.0033086]]
    
    CO2Cs3 = [[-3.3390, 1.1996, -1.0604, 0.16454],
              [0.90786, 0.086726, 0.13797, -0.035144],
              [-0.15563, -0.10292, 0.06443, -0.014128]]
    
    CO2Cs4 = [[-3.0380, 0.087994, 0.44952, -0.63679, 0.14030],
              [1.1288, -1.0822, 1.5792, -0.74749, 0.12207],
              [-0.25513, 0.045499, -0.22845, 0.16615, -0.034597],
              [0.036827, 0.040937, 0.018056, -0.031075, 0.0076346]]
    
    #we must pick which one of these that we should use, this is then used to determine M and N
    #The N and M values are determined in the emissivity function where they are called for use
    waterC = waterCs1
    CO2C   = CO2Cs4
    
    #we must now recursivly compute the natural log of the non pressure correlated emissivity values this will happen in a subfunction
    Emissivity0w = m.exp(lnEmissivity0func(waterC, tao, lamw))
    Emissivity0c = m.exp(lnEmissivity0func(CO2C, tao, lamc))

    
    
    #Now we must compute the pressure correction correlation ep/e0
    #This will be mildly inconvinient
    #we must compute everything twice once for water, and once for CO2
    PT = gas.P
    pw = PT*gas.X[gas.species_index('H2O')]
    pc = PT*gas.X[gas.species_index('CO2')]
    
    #First we must calculate an effective pressure
    PEw = PT*(1+4.8*(pw/PT)*m.sqrt(273/T))
    PEc = PT*(1+0.28*(pc/PT))

    #The function for pressure correlation is now a function of T, pL, and the new effective pressure PE
    #the maximum correction factor and position of maxima must be calculated, these are temperature dependent
    lamMw = m.log(13.2*tao**2)
    if T > 700:
        lamMc = m.log(0.225*tao**2)
    else:
        lamMc = m.log(0.054*tao**(-2))
        
    #correction factor maximum is primarily dependent on 2 functions each
    if T > 750:
        Aw = lambda t: 1.888 - 2.053*m.log(t)
    else:
        Aw = lambda t: 1.888 - 2.053*m.log(0.75)
    Bw = lambda t: 1.10*t**(-1.4)
    Ac = lambda t: 0.10*t**(-1.45)+1
    Bc = lambda t: 0.23 #spectacular that this is not even a fucntion of temperature.... what
    
    #now we must compute the maximum correction factor for both of these molecules.
    EcrMw = (Aw(tao)*PEw+Bw(tao))/(PEw+Aw(tao)+Bw(tao-1))
    EcrMc = (Ac(tao)*PEc+Bc(tao))/(PEc+Ac(tao)+Bc(tao-1))
    
    #now we need the part of this function that is dependent on the optical path length
    fpLw = m.exp(-0.5 *(lamMw-lamw)**2)
    fpLc = m.exp(-1.47*(lamMc-lamc)**2)
    
    #now we must compute the final correction factor
    Epw = ((fpLw*(EcrMw-1))+1)*Emissivity0w
    Epc = ((fpLc*(EcrMc-1))+1)*Emissivity0c
    
    
    
    
    #now that we have the pressure and temperature dependent emissivity values emperically we need to calculate the emissivity overlap from when they emit togather
    deltaE = overlap((pw/(pw+pc)), (pLc+pLw))
    
    #return the final emissivity... probably
    #return Epw + Epc - deltaE
    return [Epw, Epc, deltaE]

def lnEmissivity0func(Carr, tao, lam):
    #computes the polynomial equations needed for emissivity0
    
    #M and N values for the given coefficient arrays
    M   = len(Carr)-1
    N   = len(Carr[0])-1
    
    lnEmissivity0 = 0
    for i in range(M):
        #the first thing that we need is the a value
        a = 0
        for j in  range(N):
            a += Carr[i][j]*(tao**j) #this should compute the a coefficient
        #now this must get used to compute lnE
        lnEmissivity0 += a*lam**i
    return lnEmissivity0
    
#determines the cartesian distance between two 3 points
def cartDist(p1, p2):
    return m.sqrt((p2[0]-p1[0])**2+(p2[1]-p1[1])**2+(p2[2]-p1[2])**2)

def overlap(zeta, p):
    lamb = m.log10(p)
    return ((zeta/(10.7+101*zeta))-(0.0089*zeta**10.4))*m.log(p)**2.76









def calcBLC2():
    #Initial variable setup for important things
    BLCMdot =  IV.BLCOrificeNum * Inj.MdotSPIONLY(IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, Pc, IV.FuelTankP) #mass flow rate of boundary layer propellent
    Ul_init  = (BLCMdot/(CP.PropsSI("D", "T|liquid", IV.FuelTankT, "P", IV.FuelTankP, "Ethanol")*(IV.BLCOrificeNum*(math.pi*(IV.BLCOrificeDiameter*(in2m)/2)**2))))*(math.sin(IV.BLCAngle*d2r)) #Initial liquid coolant velocity tangent to the chamber wall
    Tv = CP.PropsSI("T", "P", ps[i], "Q", 1, "Ethanol") #Saturation temp of coolant
    Tl = IV.fuelTantT
    
    #Starting with the coolant heating methodology similar to grissom, our coolant is not at saturation temperature, so we must determine how long it takes for this to happen (long as in distance)
    x = 0 
    dx = LT/IV.CellNum
    while Tl < Tv:
        Tl = 1000
        #will need radiative and convective heat fluxes, radiative heat flux will be inconvinient, will be options to take it from logged data, and option to compute it wil ray tracing, but will use simpler correlations.abs
        
        
        
        x = x + dx #get the next x value for distance downstream
    
    

    return 0

#Initial variabl setup for some important starting values

CombustionGas = equalibrateHPWrapper(IV.Pc, IV.OF, 10**(-3))
#print(CombustionGas.X)
#print(getRay(0, (75*(m.pi/180)), (0*(m.pi/180))))
#print(getEmittanceAtRay(0, (75*(m.pi/180)), (0*(m.pi/180)), ceaout[0]))
print(getEmittanceAtRay(0, (0*(m.pi/180)), (0*(m.pi/180)), CombustionGas))


#Test Code down here, to demonstrate the validity of the sub functions, primarily the ones
#that have graphs to compare against. Specifically as listed below
#Emissivity ratio
#Emissivity ratio maximum
#Emissivities of water vapor
#Emissivities of carbon dioxide
#Overlap delta

if IV.EmittanceTest:
    #first we are going 
    print(1)

