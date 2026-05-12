# Right now this is incredibally simplue, just uses the SPI model, wont work for all props work for now though
"""
@author: phineas
"""

import CoolProp.CoolProp as CP
import math

in2meter = 0.0254
PSI2PA = 6894.76


def MdotSPIONLY(Cd, D, prop, TProp, PChamber, PTank):
    A = math.pi * 0.25 * ((D * in2meter) ** 2)
    rho0 = CP.PropsSI("D", "T", TProp, "Q", 0, prop)
    Mdot = Cd * A * ((2 * rho0 * ((PTank * PSI2PA) - PChamber)) ** 0.5)
    return Mdot


def Mdot(Cd, D, prop, TProp, PChamber, PTank):  # Currently Nonfunctional
    Pvap = CP.PropsSI("P", "T", TProp, "Q", 1, prop)
    A = math.pi * 0.25 * ((D * in2meter) ** 2)
    rho0 = CP.PropsSI("D", "T", TProp, "Q", 0, prop)

    h0 = CP.PropsSI("H", "T", TProp, "Q", 0, prop)
    s0 = CP.PropsSI("S", "T", TProp, "Q", 0, prop)
    rho1 = CP.PropsSI("D", "S", s0, "P", PChamber, prop)
    h1 = CP.PropsSI("H", "S", s0, "P", PChamber, prop)

    k = ((PTank - PChamber) / (Pvap - PChamber)) ** 0.5

    MdotHEM = Cd * A * rho1 * ((2 * (h0 - h1)) ** 0.5)
    MdotSPI = Cd * A * ((2 * rho0 * ((PTank * PSI2PA) - PChamber)) ** 0.5)
    print("PTank" + str(PTank))

    if PChamber > Pvap:
        print("Pc = " + str(PChamber) + " : " + "Pvap = " + str(Pvap))
        return (k / (1 + k)) * MdotSPI + (1 / (1 + k)) * MdotHEM
    else:
        print("here2")
        return MdotSPI
