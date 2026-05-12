import math
import numpy as np
import cantera as ct
import CoolProp.CoolProp as CP
from .Utils.cea import runCEA, AxialValues
from .Utils.EngineGeometry import RatL, LT
from .Inputs import InputValues as IV
from .Utils import Injector as Inj
#from .Utils.Meshing import MeshGen as MG



AxialDistances = np.linspace(0, LT, IV.CellNum)
RadiusVal = [0.0]*IV.CellNum
#Creates radial values for each axial length value
for i in range(IV.CellNum):
    RadiusVal[i] = RatL(AxialDistances[i])

#for i in range(IV.CellNum):
    #print(AxialDistances[i])


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
print("F : " + str(v*(FuelMdot+OxMdot)*0.9+((math.pi*(RadiusVal[IV.CellNum-1]*0.0254)**2)*(val[1][IV.CellNum-1]-101000))))
print("V : " + str(v))



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
