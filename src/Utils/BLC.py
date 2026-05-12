# No like truly this is a lot
# fmt: off
import InputValues as IV
from cea import runCEA, AxialValues
import math as m
import Injector as Inj
from EngineGeometry import LT, Dc, RatL
import mathfunctions as mf
import cantera as ct
import CoolProp.CoolProp as CP
from Bisect import Bisect


pa2atm = 9.86923 * 10 ** (-6)
in2m = 0.0254

#This function calculates the emittance of the gas for some point in the rocket engine based on the diameter and the CombusionGas
def ε(CombusionGas, R, x):
    # effective length through gas
    LeffD = 0.95*2*R
    LeffU = 0.95*2*((4*x)/(2*R+4*x))
    Leff = 0.5*(LeffD+LeffU)*IV.Aw**(-0.85)

    # Partial Pressures calculated form mole fraction and total Pressures
    χC = CombusionGas.mole_fraction_dict()["CO2"]
    χH = CombusionGas.mole_fraction_dict()["H2O"]
    ppC = χC * CombusionGas.P * pa2atm
    ppH = χH * CombusionGas.P * pa2atm

    # Optical Densities Based on the effective length through the gas and the partial pressure, should probably be a path integral
    ρOptC = ppC * Leff
    ρOptH = ppH * Leff

    T = CombusionGas.T
    # These are the 4 tempature based parabolas from the table in Grisson for the c and n coefficients this should probably be replaced with a general polyfit function
    cC = (2.5 * 10 ** (-8)) * T**2 + (-0.00005) * T + (0.075)
    nC = (0) * T**2 + (0) * T + (0.6)
    cH = (2.075 * 10 ** (-7)) * T**2 + (0.0001125) * T + (-0.155)
    nH = (-1.2 * 10 ** (-7)) * T**2 + (0.00056) * T + (0.01)

    # Calculates uncorrecected emittance values without the Kp pressure correction
    εfC = 0.231
    εfH = 0.825
    εC = εfC * (1 + (ρOptC / cC) ** (-nC)) ** (-1 / nC)
    εH = εfH * (1 + (ρOptH / cH) ** (-nH)) ** (-1 / nH)

    # Pressure correction factors for the given emittance values at 1 atm
    KpC = 10 ** ( 0.036 * ρOptC ** (-0.489) * (1 + (2 * m.log(CombusionGas.P*pa2atm, 10)) ** (-100 * ppC)) ** (-1 / (100 * ppC)))

    c1 = 0.26 + 0.74 * m.exp(-2.5 * ρOptH)
    c2 = 0.75 + 0.31 * m.exp(-10 * ρOptH)
    KpH = 1 + c1 * ( 1 - (m.exp((1 - CombusionGas.P*pa2atm * (1 + χH)) / c2)))

    # Calculates the emittance correction factor to handle ovlerlap in the 2 specturms
    n = 5.5 * (1 + (1.09 * (ρOptC + ρOptH)) ** (-3.88)) ** (-1 / 3.88)
    Kx = 1 - (abs((2 * χH)/( χH + χC) - 1)) ** n
    Δε = 0.0551*Kx*(1-m.exp(-4*(ρOptC+ρOptH)))*(1-m.exp(-12.5*(ρOptH+ρOptH)))

    # Return the final emittance (Probably)
    return (KpC*εC) + (KpH*εH) - (Δε)

#Funciton to calculate Prandtl number of a cantera gas
def getPrandtl(gas):
    return (gas.cp*gas.viscosity)/gas.thermal_conductivity

#Start of BLC calc finction that iterates down the chamber length
def calcBLC():
    #Starting Paramaters
    x = 0
    dL = LT/IV.CellNum
    Lold = 0
    Rold = (Dc*in2m)/2
    Tc = IV.FuelTankT
    Tca = [0.0]*IV.CellNum

    #Combustion Paramaters
    ceaOut = runCEA()
    CombustionGas = ceaOut[0]
    γ = CombustionGas.cp/CombustionGas.cv
    cps = CombustionGas.cp

    #Axial Values along engine
    Values = AxialValues(CombustionGas.T, CombustionGas.P, CombustionGas.density, CombustionGas)
    Ts = Values[0]
    ps = Values[1]
    ρs = Values[2]
    Tr = Values[3]
    Ms = Values[4]

    #Calculates Velocity Based on engine data in values and Mach number also does radii and lengthto save a loop
    L = mf.linspace(0, LT, IV.CellNum)
    Rad = [0.0]*IV.CellNum
    Us = [0.0]*IV.CellNum
    for i in range(IV.CellNum):
        Rad[i] = RatL(L[i])*in2m
        L[i] *= in2m
        Us[i] = Ms[i]*(γ*Ts[i]*(ct.gas_constant/CombustionGas.mean_molecular_weight))**0.5

    #This is layed out this way in preperation for multiphase injection
    BLCMdot = IV.BLCOrificeNum * Inj.MdotSPIONLY( IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, ceaOut[0].P, IV.FuelTankP)
    BLCMdotL = BLCMdot

    print("BLC: ---------------- " + BLCMdot)
    #Boltzman const
    σ = 5.67*(10**(-8))

    #Coolant Flow length per circumfrence
    Γ = BLCMdot/(Rad[0]*2*m.pi)

    for i in range(IV.CellNum):
        #Calculates Contour length and axial length
        dx = (((L[i]-Lold)**2)+(Rad[i]-Rold)**2)**0.5
        Lold = L[i]
        Rold = Rad[i]
        x += dx

        #Molecular weights of stream gas and Coolant
        Mwc = CP.PropsSI("M", "P", ps[0], "T", IV.FuelTankT, "Ethanol")
        Mws = CombustionGas.mean_molecular_weight/1000

        #Checks to see if BLC liquid still remains
        if Γ > 0:
            #This is an implicit formulation derrived from Grissom, Hopefully I can work this out for someone else:
            #This was pain and this still might eb the wrong function
            Tv = CP.PropsSI("T", "P", ps[i], "Q", 1, "Ethanol") #Saturation temp
            xe = 3.53*(Rad[i]*2)*(1+(x/(3.53*Rad[i])+(10**(-5)))**(-1.2))**(-1/1.2)
            Gch = ρs[i]*Us[i]
            Tm = 0.5*(Ts[i]+Tv)

            #Actuall function def for bisecting
            def Ulf(Ul):
                const = ((0.0592)/(2*0.023))**5
                a = (Gch*(Ts[i]/Tm)*((Us[i]-Ul)/Us[i]))**2
                b = const*(1/(xe*Rad[i]*2))
                return a - b
            Ul = Bisect(Ulf, 0, Us[i], 10*(-20))

            #This if statement makes a better starting Ul value
            if i == 0:
                Ul = (BLCMdot/IV.BLCOrificeNum)/(CP.PropsSI("D", "T|liquid", IV.FuelTankT, "P", IV.FuelTankP, "Ethanol")*(m.pi*((IV.BLCOrificeDiameter*in2m)/2)**2))
                print("U:" + str(Ul))

            #Ul = 10
            print(Ul)

            #This section now goes to calculate gas -> liquid heat transfer coefficient
            Kt = 1+4*IV.et #Turbulence Correction factor
            ug = CombustionGas.viscosity
            G = Gch*(Ts[i]/Tm)*((Us[i]-Ul)/Us[i]) 
            Rex = G*(xe/ug) + 10**(-5) #Renolds number basec on effective contour length
            Cf = 0.0592*Rex**(-0.2) #Skin friction coefficiet
            Pr = getPrandtl(CombustionGas)
            St0 = 0.5*Cf*Pr**(-0.6) #Stanton Number
            h0 = Kt*G*St0*CombustionGas.cp

            #If the coolent is not at saturation temp coolant temp goes up
            if Tc < Tv:
                #Transperation is ignored, there is none happening yet
                Qconv = h0*(Tr[i]-Tc)
                Qrad = σ*IV.Aw*ε(CombustionGas, Rad[i], x)*((Ts[i]**4)-(Tv**4))
                Qtot = Qconv+Qrad

                #Finds coolant specific heat capacity and calcs temp rise per dx
                Coolcp = CP.PropsSI("C", "T", Tc, "P", ps[i], "Ethanol")
                dT = (Qtot/(Γ*Coolcp))*dx
                Tc += dT

            #Else, coolant is at saturation temp, so we need to eavaporate some now
            else: 
                #Calculates latent heat of vaporization
                LHVc = CP.PropsSI("H", "P", ps[i], "Q", 1, "Ethanol")-CP.PropsSI("H", "P", ps[i], "Q", 0, "Ethanol")

                #Calculate stream -> Liquid heat transfer coefficient now with Transperation
                if Mwc < Mws:
                    a = 0.6
                else:
                    a = 0.3

                #Function defenition for H
                def H(h):
                    Qconv = h*(Tr[i]-Tc)
                    Qrad = σ*IV.Aw*ε(CombustionGas, Rad[i], x)*((Ts[i]**4)-(Tv**4))
                    Qtot = Qconv+Qrad

                    Km = (Mws/Mwc)**a
                    return (cps*Km*(Qtot/LHVc))/h

                #Defenition for gas heat transfer coeff with transperation
                def hf(h):
                    return (h0*m.log(1+H(h)))/H(h) - h

                #Uses the function defined to get h
                h = Bisect(hf, 10**(-20), 10**10, 10**(-20))

                Qconv = h*(Tr[i]-Tc)
                Qrad = σ*IV.Aw*ε(CombustionGas, Rad[i], x)*((Ts[i]**4)-(Tv**4))
                Qtot = Qconv+Qrad

                #Mdot vapor, this is mass increase of vapor fraction and decrease in liquid fraction
                ṁvap =  Qtot/LHVc
                dΓ = -ṁvap*dx #Evaporation rate per contour length
                Γ += dΓ

                #Uh uh we ran out of film Coolant
                if Γ <= 0:
                    Γ = 0

                    #Calculates Initial entrained gas flow for gas cooling section
                    G_local = ρs[i]*Us[i]
                    mus = CombustionGas.viscosity
                    K = G_local*(mus**0.25)*((BLCMdot/(Rad[i]*2*m.pi))**(-1.25))
                    Xi = K*x
                    MC_bl

        #Its gas time now
        else:
            print("length " + str(L[i]/in2m))
        Tca[i] = Tc
                 
    return Tca

ceaOut = runCEA()
BLCMdot = IV.BLCOrificeNum * Inj.MdotSPIONLY( IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, ceaOut[0].P, IV.FuelTankP)
CombustionGas = ceaOut[0]
#print(emittance(CombustionGas, 3*0.0254))
#print(((calcBLC()[249]*ceaOut[1])+(BLCMdot*(0.8*calcBLC()[249])*0.9)/9.81)/(BLCMdot+ceaOut[1]))
#print(RatL(LT)/12)
#Values = AxialValues(CombustionGas.T, CombustionGas.P, CombustionGas.density, CombustionGas)
#print(Values[1])
#print(CombustionGas.report())
#print(calcBLC())
