# CEA calculations for rocket engines.
# fmt: off
"""
@author: phineas
"""
#External modules
import math as m
import numpy as np
import CoolProp.CoolProp as CP
import cantera as ct

#Internal modules
from .EngineGeometry import RatL, Dt, LT, Lt
from ..Inputs import InputValues as IV
from ..Utils import Injector as Inj
from .Bisect import Bisect
from ..Helpers.PropName import conCool, conCant

in2m = 0.0254 #inch 2 meter conversion, should be moved to a different file
psi2pa = 6894.71 #Psi to pascals

def runCEA():
    #Starting conditions for iteration
    Pc = IV.AmbP 
    Tc = IV.AmbT
    Tol = 10 ** (-4)
    RelError = 1
    Damp = 0.6

    #Creates CombustionGas
    CombustionGas = ct.Solution(IV.yaml)
    #print(CombustionGas.report())

    i = 0
    while (RelError > Tol) & (i < 250):
        i += 1  # This is here to ensure no infinite loops
        PcOld = Pc #Keeps track of old Pc to find reletive error
        TcOld = Tc #Keeps track of old Tc to find reletive error


        # Gets Mdot for Both fuel and Ox sides all orifices, also total Mdot,
        FuelMdot = IV.FuelOrificeNum * Inj.MdotSPIONLY( IV.FuelOrificeCd, IV.FuelOrificeDiameter, IV.Fuel, IV.FuelTankT, Pc, IV.FuelTankP)
        BLCMdot = IV.BLCOrificeNum * Inj.MdotSPIONLY( IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, Pc, IV.FuelTankP)
        OxMdot = IV.OxOrificeNum * Inj.MdotSPIONLY(IV.OxOrificeCd, IV.OxOrificeDiameter, IV.Ox, IV.OxTankT, Pc, IV.OxTankP)
        Mdot = FuelMdot+OxMdot+BLCMdot
        #print(FuelMdot)
        #print(OxMdot)
        
        #Calculates OF ratio, technically not efficient to have it here or caculated this way, but eh.
        OF = OxMdot / FuelMdot

        #get the combustion gases at the current chamber pressure and OF ratio
        CombustionGas = equalibrateHPWrapper(Pc, OF, 10**(-4))
        
        #Calculates new ChamberPressure from the new tempature info, and gas propertys
        γ = CombustionGas.cp/CombustionGas.cv
        R = ct.gas_constant/CombustionGas.mean_molecular_weight
        Tc = CombustionGas.T
        Pc = ChamberPressure(Tc, Mdot, γ, R)
        #print((γ*R*Tc)**0.5)
 
        #Calculates Max reletive error between ChamberPressure and ChamberTempature
        RelError = max((abs(Pc-PcOld))/(Pc),(abs(Tc-TcOld))/(Tc))

        #Calculates the chamge in ChamberPressure and ChamberTempature, makes sure its not so large it just overshoots everything and the engine "explodes"
        Pc -= Damp*(Pc-PcOld)
        Tc -= Damp*(Tc-TcOld)
        #print(Mdot)
        #print("OxMdot: " + str(OxMdot+FuelMdot+BLCMdot))
        #print("FuelMdot: " + str(FuelMdot))
    #print(CombustionGas.P)
    #print(OF)
    return CombustionGas, Mdot



#Function for determining chamber pressure
def ChamberPressure(Tc, mdot, gamma, R):
    exp = -0.5 * ((gamma + 1) / (gamma - 1))
    u = (gamma * R * Tc) ** 0.5
    At = m.pi * ((Dt * in2m / 2) ** 2)
    MdotTerm = mdot / (At * gamma)

    Pc = MdotTerm * u * (2 / (gamma + 1)) ** exp
    return Pc

#Cantera does not support fluids so, we must input the correct enthalpy for our cold fluids to have the correct energy for combustion
#However in many cases with fluids, particularly cold fluids the total enthalpy of the ideal gas would have a negetive temperature, cantera does not like this
#So we must find a workaround, here we assume that an HP equiliberation would somewhat believably keep constant enthalpy and pressure
#We use TP equiliberation to find some temperature that has the expected enthalpy, its less efficient but we never have to set a gas with low chemical potential to
#Our inlet enthalpy this method seems to agree, for the nice cases, with directly setting enthalpy getting a 0.8K gas and equalibrating
def equalibrateHPWrapper(Pc, OF, tol):
    #First we must compute the input enthalpy from our propellents... this is odd as cantera uses a vastly different enthalpy system than coolprop but its what we must do
    #Cantera uses a reference zero for its species at 1atm and 0c so we rezero out coolprop values to this
    HF = CP.PropsSI("H", "P", Pc, "T", IV.FuelTankT, conCool(IV.Fuel))-CP.PropsSI("H", "P", 101335, "T", 273.15, conCool(IV.Fuel))
    Ho = CP.PropsSI("H", "P", Pc, "T", IV.OxTankT, conCool(IV.Ox))-CP.PropsSI("H", "P", 101335, "T", 273.15, conCool(IV.Ox))
    
    #We then find the enthalp of the propellent mix
    Hin =  (HF * (1 / (1 + OF))) + (Ho * (OF / (1 + OF)))
    
    #Unfortunatly cantera uses absolute enthalpy so, even if we are at zero sensible enthalpy, there will still be heat of formation, so we must find that part too
    TestGas = ct.Solution(IV.yaml)
    #print(conCant(IV.Fuel, IV.Ox, OF))
    TestGas.TPY = 273.15, 101325, conCant(IV.Fuel, IV.Ox, OF)
    
    #We can add them to get the total needed enthalpy of the gas
    EnthalpyTarget = TestGas.h + Hin

    #If we are assuming that we want to do HP constant combustion here, we would need to set this in a gas, however that may not work
    #So what we will do is try varous temperatures untill one of them gives us the desired enthalpy with a bisect method on a function shown below
    Treaction = Bisect(equalibrateHP, 10, 5000, tol, EnthalpyTarget, Pc, OF)
    
    #now that we have the known temperature that produces our desired enthalpy we can construct our gas
    CombustionGas = ct.Solution(IV.yaml)
    CombustionGas.TPY = Treaction, Pc, conCant(IV.Fuel, IV.Ox, OF)
    CombustionGas.equilibrate("TP")

    return CombustionGas



#This function makes a gas at some pressure and temperature, and returnes the equalibrium enthalp for the fuel mixture defined
def equalibrateHP(T, H, Pc, OF):
    #First we must create the gas at temperature pressure and of
    CombustionGas = ct.Solution(IV.yaml)
    CombustionGas.TPY = T, Pc, conCant(IV.Fuel, IV.Ox, OF)
    
    #We now find the equilibrium and get the enthalp
    CombustionGas.equilibrate("TP")
    return (H-CombustionGas.h)



#this function takes in baisc charictaristics and outputs engine geometris and such
def basicEngineChar(Pc, Pe, OF, F):
    #Produces the needed combustion gas for finding out engine charictaristics
    CombustionGas = equalibrateHPWrapper(Pc, OF, 10**(-5))
    
    #We we get the 
    y = CombustionGas.cp/CombustionGas.cv
    R = ct.gas_constant/CombustionGas.mean_molecular_weight
    Tc = CombustionGas.T
    
    #Ideal characteristic velocity
    exp = (y+1)/(y-1)
    cStar = ((m.sqrt(y*R*Tc))/(y*m.sqrt((2/(y+1))**exp)))*IV.cStarEff
    
    #Thrust Coeff: consider, the nozzle is iedally expanded, so  p2-p3=0 in sutton:3-30
    frac1 = (2*y**2)/(y-1)
    frac2 = (2/(y+1))
    frac3 = (y+1)/(y-1)
    CF = m.sqrt(frac1*(frac2**frac3)*(1-(101325/Pc)**((y-1)/y)))
    
    #Compute mass flow from thrust and exaust velocity we assume that it is ideally expanded so thrust is Ue*mdot
    Ue = (cStar*CF)
    mdot = F/Ue #Mass flow rate in Kg/s
    
    #we can also find the area of the throat from the cmaber pressure and mass flow sutton 3-24
    frac1 = 2/(y+1)
    frac2 = (y+1)/(y-1)
    fracRHS = (m.sqrt(frac1**frac2))/(m.sqrt(y*R*Tc))
    RHS = fracRHS*y*Pc
    At = mdot/RHS #Area of the throat
    
    #now we should also get exis area ratio yipeee but first we need the mach number in the chamber
    #we will use that to get stagnation values
    Ac = m.pi*(IV.Dc*in2m/2)**2 #Chamber area in square meters
    Uc = mdot/(Ac*CombustionGas.density)
    Mc = Uc/m.sqrt(y*R*Tc)
    P0 = Pc/((1+((y+1)/2)*Mc**2)**(-y/(y-1))) #Stagnation pressure
    
    #We need mach number at the exit, so we must determine the temperature at exit
    T0 = Tc/((1+((y+1)/2)*Mc**2)**(-1))
    
    #Use T0 to find exit temperature from pressure ratio knowing that exit pressure is 1atm
    Te = T0*(Pe/Pc)**((y-1)/y)
    Me = Ue/(m.sqrt(y*R*Te))
    
    #now we calculate expansion ratio
    frac1 = 1+(((y-1)/2)*(Me**2))
    frac2 = 1+(((y-1)/2))
    frac3 = (y+1)/(y-1)
    ExpR  = (1/Me)*m.sqrt((frac1/frac2)**frac3)
    
    #Calculate exit area
    Ae = At*ExpR
    
    #Return All of the desirec rocket engine values in a list
    return [cStar, CF, mdot, At, ExpR, Ae]



#This function fionds the axial values for several things, temp pressure adiabatic wall temp, etc
def AxialValues(Tc, pc, ρc, CombustionGas):
    #Crates array of Lengths and an array of radisus
    AxialDistances = np.linspace(0, LT, IV.CellNum)
    RadiusVal = [0.0]*IV.CellNum

    #Creates radial values for each axial length value
    for i in range(IV.CellNum):
        RadiusVal[i] = RatL(AxialDistances[i])
    
    #Creates the arrays for Engine values along the axial length
    Ms = [0.0]*IV.CellNum #Mach number in stream
    Ts = [0.0]*IV.CellNum #Temp in stream
    Tr = [0.0]*IV.CellNum #Recovery temp at "wall"
    ps = [0.0]*IV.CellNum #Pressure in stream
    ρs = [0.0]*IV.CellNum #Density in stream


    #Gamma from CombustionGas
    γ = CombustionGas.cp/CombustionGas.cv

    #Loops through all AxialDistances to get the values there
    for i in range(IV.CellNum):
        #The Mach number function here gets redefined at every iteration, this should be moved and the Mach number function should
        #just need 2 input values (x, A) but I have not looked into variable passthrough yet
        A = m.pi*RadiusVal[i]**2
        def MachNumber(x):
            At = 0.25*m.pi*IV.Dt**2
            exp = (γ+1)/(γ-1)
            mult = 0.5*(γ-1)

            #Returnes the mach number at the given area ratio here
            return ((1/x)*(((1+(mult*x**2))/(1+mult))**exp)**0.5)-(A/At)
    
        if AxialDistances[i] < Lt:
            Ms[i] = Bisect(MachNumber, 0.001, 0.99999, 10**-20)
        else:
            Ms[i] = Bisect(MachNumber, 1.00001, 5, 10**-20)
        
        #Uses Isentropic FLow equationjs to calculatie Stream tempature, Pressure, and density
        Ts[i] = Tc*(1+(0.5*(γ-1))*Ms[i]**2)**(-1)
        ps[i] = pc*(1+(0.5*(γ-1))*Ms[i]**2)**((-γ)/(γ-1))
        ρs[i] = ρc*(1+(0.5*(γ-1))*Ms[i]**2)**((-1)/(γ-1))

        #Calculates the recovery tempature
        Pr = (CombustionGas.cp*CombustionGas.viscosity)/(CombustionGas.thermal_conductivity) 
        r = Pr**0.33
        mult = (0.5*(γ-1))*Ms[i]**2
        Tr[i] = Tc*((1+r*(mult))/(1+mult))

    return Ts, ps, ρs, Tr, Ms