# I think this will end up being a quazi 1D CFD for providing more precice chemical reactions along the nozzle
# from EngineGeometry import RatL, LT, Lt
import InputValues as IV

# from ChamberPressure import ChamberPressure
##import Injector as Inj
# import mathfunctions as mf
import cantera as ct

# import CoolProp.CoolProp as CP
# import math
# from Bisect import Bisect

# The cantera objects will only store reaction data and species ratios, mdot will have to be used
# to determine species, hopefully this does not fuck with species transport like I think it might


# Simulates a single timestep of the gas, turns out you have to make a ReactorNet
# of size 1 to get finite rate chemistry to happen
def timestepGas(CombustionGas, t, dt):
    # Creats a Cantera reactor object
    react = ct.IdealGasReactor(CombustionGas)

    # Creats a reactor network to simulate time its 1 reactor in the network, thank you cantera
    reactNet = ct.ReactorNet([react])

    # For some unknows reason to me Cantera reactors DONT LIKE TIMESTEPS, so this has to sim to time+timestep
    reactNet.advance(t + dt)
    return CombustionGas


OF = 1.6
TRef = 3000
Pc = 100000
CombustionGas = ct.Solution(IV.yaml)
CombustionGas.TPY = TRef, Pc, "O2:" + str(OF) + ", C2H5OH:1"
# print(CombustionGas.forward_rates_of_progress)
# print(CombustionGas.chemical_potentials)

# print(timestepGas(CombustionGas, 0, 0.000001).report())
# CombustionGas.equilibrate("TP")
print("Species: " + str(CombustionGas.n_species))

t = 0
dt = 10 ** (-10)
tend = 10 ** (-9)
T = [0.0] * int(tend / dt)
for i in range(int(tend / dt)):
    T[i] = CombustionGas.T
    CombustionGas = timestepGas(CombustionGas, t, dt)
    t += dt
    if i % 1000 == 0:
        print(i)

for i in range(len(T)):
    print(T[i])
