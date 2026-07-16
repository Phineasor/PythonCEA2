import math as m
import numpy as np
import cantera as ct
import CoolProp.CoolProp as CP

from .Utils.cea import runCEA, AxialValues, basicEngineChar
from .Utils.EngineGeometry import RatL, LT
from .Inputs import InputValues as IV
from .Utils import Injector as Inj
from .Helpers import GraphGen
from .Helpers.PropName import conCool, conCant
from .Utils.BLC2 import calcBLC2 #Temporarily Disabled
#from .Utils.Meshing import MeshGen as MG

in2m = 0.0254
psi2pa = 6894.71

AxialDistances = np.linspace(0, LT, IV.CellNum)
RadiusVal = [0.0]*IV.CellNum
#Creates radial values for each axial length value
for i in range(IV.CellNum):
    RadiusVal[i] = RatL(AxialDistances[i])

#for i in range(IV.CellNum):
    #print(AxialDistances[i])



if IV.Physical:
    print("--------------------------------------------------------------------------------")
    print("Physical driven rocket engine values")
    ceaOut = runCEA()
    γ = ceaOut[0].cp/ceaOut[0].cv
    R = ct.gas_constant/ceaOut[0].mean_molecular_weight
    val = AxialValues(ceaOut[0].T, ceaOut[0].P, ceaOut[0].density, ceaOut[0])
    print("ExitPressure: " + str(val[1][IV.CellNum-1]))
    print("ExitTemp: " + str(val[0][IV.CellNum-1]))

    v = ((γ*R*val[0][-1])**0.5*val[4][-1])
    BLCMdot = IV.BLCOrificeNum * Inj.MdotSPIONLY( IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, ceaOut[0].P, IV.FuelTankP)
    FuelMdot = IV.FuelOrificeNum * Inj.MdotSPIONLY( IV.FuelOrificeCd, IV.FuelOrificeDiameter, IV.Fuel, IV.FuelTankT, ceaOut[0].P, IV.FuelTankP)
    OxMdot = IV.OxOrificeNum * Inj.MdotSPIONLY(IV.OxOrificeCd, IV.OxOrificeDiameter, IV.Ox, IV.OxTankT, ceaOut[0].P, IV.OxTankP)
    print("OF: " + str(OxMdot/FuelMdot))
    print("BLCMdot: " + str(BLCMdot))
    print("BLCMdot/BLCMdot+F: " + str(BLCMdot/(FuelMdot+BLCMdot)))
    print("pc: " + str(ceaOut[0].P*0.000145038))
    print("Tc: " + str(ceaOut[0].T))
    print("Te: " + str(val[0][-1]))
    print("F : " + str(v*(FuelMdot+OxMdot)*0.9+((m.pi*(RadiusVal[IV.CellNum-1]*0.0254)**2)*(val[1][IV.CellNum-1]-101000))))
    print("Ue: " + str(v))

print("--------------------------------------------------------------------------------")
print("Standard Chamber Size, Pressure and Thrust driven engine paramaters")

List = basicEngineChar(IV.Pc, IV.AmbP, IV.OF, IV.F) #gets all of the known values and uses them to produce the needed engine params
print("cStar: " + str(List[0]))
print("Thrust Coefficient: " + str(List[1]))
print("Total Mass flow rate [Kg/s]: " + str(List[2]))
MdotF = (1 / (1 + IV.OF))*List[2]
MdotO = (IV.OF / (1 + IV.OF))*List[2]
print("Fuel Mass flow [Kg/s]: " + str(MdotF))
print("Ox Mass flow rate [Kg/s]: " + str(MdotO))
print("Throat Area [m^2]: " + str(List[3]))
print("Throat diameter [in]: " + str(m.sqrt(List[3]/m.pi)*2/in2m))
print("ExpRatio: " + str(List[4]))
print("Exit diameter [in]: " + str(m.sqrt(List[5]/m.pi)*2/in2m))
Vc = (List[3]/(in2m**2))*IV.LStar
print("Chamber Length [in]: " + str(Vc/(m.pi*(IV.Dc/2)**2)))

print(".")
#need oriface area 
sqrtF = m.sqrt(2*CP.PropsSI("D", "T", IV.FuelTankT, "P", IV.Pc, conCool(IV.Fuel))*(IV.FuelTankP*psi2pa-IV.Pc))
sqrtO = m.sqrt(2*CP.PropsSI("D", "T", IV.OxTankT, "P", IV.Pc, conCool(IV.Ox))*(IV.OxTankP*psi2pa-IV.Pc))
AF = MdotF/(0.7*sqrtF)
AO = MdotO/(0.7*sqrtO)
print("Area of fuel oriface [in^2]: " + str(AF/(in2m**2)))
print("Area of ox   oriface [in^2]: " + str(AO/(in2m**2)))
print(CP.PropsSI("D", "T", IV.FuelTankT, "P", IV.Pc, conCool(IV.Fuel)))
print(CP.PropsSI("D", "T", IV.OxTankT, "P", IV.Pc, conCool(IV.Ox)))


#print(CP.PropsSI("T", "P", CombustionGas.P, "Q", 1, "Ethanol"))
#print(CP.PropsSI("T", "P", CombustionGas.P, "Q", 0.5, "Ethanol"))
#print(ceaOut[0]())

#print(AxialValues(CombustionGas.T, CombustionGas.P, CombustionGas.density, CombustionGas)[3])
#print(ceaOut[0].report())
#print(str(Pc/Inj.PSI2PA)+" : "+str(Pc))
#print("Pr = "+str(Pr))
#print("mu = "+str(CombustionGas.viscosity))
#print("Cp = "+str(CombustionGas.cp))
#print(CombustionGas.thermal_conductivity)
#print("Ox = "+str(264.172*10*OxMdot/CP.PropsSI("D", "T|liquid", IV.OxTankT, "P", IV.OxTankP*psi2pa, "O2")))
#print("Fuel = "+str(264.172*10*(FuelMdot+BLCMdot)/CP.PropsSI("D", "T|liquid", IV.FuelTankT, "P", IV.FuelTankP*psi2pa, "Ethanol")))
#print("Fuel gal/s = "+str(264.172*(FuelMdot+BLCMdot)/CP.PropsSI("D", "T|liquid", IV.FuelTankT, "P", IV.FuelTankP*psi2pa, "Ethanol")))
