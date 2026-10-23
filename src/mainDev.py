from pathlib import Path
import matplotlib.pyplot as plt
import cantera as ct
from .Inputs import InputValues as IV
from .Helpers.PropName import conCool, conCant
import numpy as np
import math as m
from .Utils.EngineGeometry import RatL, Dt, LT, Lt
from .Utils.Bisect import Bisect, Derivative
from scipy.integrate import ode
from scipy.integrate import solve_ivp
from joblib import parallel_config, Parallel, delayed
from .Utils.cea import runCEA, AxialValues, basicEngineChar, equalibrateHPWrapper
from .Utils import Injector as Inj



#   $$E = \frac{\sum_{\mu\nu} P_{\mu\nu} (h_{\mu\nu} + F_{\mu\nu})}{2}$$


# We will start with a "Simple" Example removing reacting features, and making it a one species gas of arbitrary property and no reaction rate
# z is our singular independent variable

# our state vector for this run will be y = [x, u, T, rho, p, M, YN] Eventually will need all
# We will need to define the functions for the subonic and  supersonic conditions

# its going to be much better to define all of our stiff G functions first, and a function for the simmation term...a


# Subsonic function
def funSubsonic(z, y):
    #unpack the state vector
    x = y[0]
    u = np.exp(y[1])
    T = np.exp(y[2])
    rho = np.exp(y[3])
    p = np.exp(y[4])
    M = y[5]
    YN = y[6:]  # These are mass fraction



    # Decault gas for geeding to funcitons when one is needed
    DgasTPY = ct.Solution(IV.yaml)

    # We will need the vector later multiple times, so it gets computed here once
    #vector = makeVector(YN, DgasTPY)
    vector = YN

    # we need to calculate some standard gas properties from our system
    Wi = np.array(DgasTPY.molecular_weights)


    # mean molecular weight aswell
    #print([T, p, vector])
    DgasTPY.TPY = T, p, vector
    W = DgasTPY.mean_molecular_weight
    #print(W)


    # for reaction rates we will need a gas at the current temperature and pressure
    DgasTPY.TPY = T, p, vector
    widot = np.array(DgasTPY.net_production_rates)


    # current specific gas contant
    R = ct.gas_constant / W

    # specific heats at constant pressure, at temp and pressure, and ratio of specific heats
    cpi = np.array(DgasTPY.partial_molar_cp) / Wi
    cp = DgasTPY.cp_mass
    cv = DgasTPY.cv_mass
    gamma = cp / cv

    # enthalpies too
    hi = np.array(DgasTPY.partial_molar_enthalpies) / Wi

    # other length dependent funcitons we should only call once
    A = Af(x)
    qw = qwf(x)
    Pw = Pwf(x)
    dAdx = Derivative(x, Af, 10**(-3),)
    Cf = Cff(x)
    Dh = Dhf(x)

    # we need mdot aswell
    mdot = rho * A * u
    
    # The first thing that we need to do is compute the gas properties from species data
    # To start this we need the rate of change for every species per distance
    dYidx = dYidxf(u, rho, Wi, widot)

    # rest of the mixture props
    dWdx = dWdxf(YN, dYidx, Wi, W)
    dRdx = dRdxf(dWdx, W, R)
    
    
    
    #there is deviation from the primary source here, Cha states that this is a caloricaly perfect gas, however cantera no do that cp is temperature dependent
    #This means that we need to redetermine dcpdx assuming that cp is a function of temperature so we must do d/dx ourselves
    
    #because this is a reacting flow, enthalpy is not constant, so we need to compute the change in enthalpy over distance
    dh0dx = dh0dxf(mdot, qw, Pw)
    
    #this is also dependent on ho and sumterm so those go ehre too
    sumterm = sumtermf(hi, dh0dx, dYidx)

    #we need dTdx for this formulation
    GT = GTf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx)
    dTdx = T*(GT/(1-M**2))
    
    #we also need dcp/dT this needs to be found numerically
    dcpdT = Derivative(T, cpTf, 10**(-3),  p, vector, DgasTPY)
    dcpdx = dcpdxf(dYidx, cpi, dcpdT, dTdx)
    

    # ratio of specific heats type shit
    dgammadx = dgammadxf(T, p, gamma, cp, R, dcpdx, dRdx)

    # we now need to determine the G functions
    Gu = Guf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx)
    #GT calculated earlier
    Grho = Grhof(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx)
    Gp = Gpf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx)
    GM = GMf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx, dgammadx)

    # [x, u, T, rho, p, M, YN]
    # now we need to compute the derivatives of all our state variables with respect to z
    dxdz = abs(1 - M**2)
    
    #These variables were transformed into log form so that the solver is more stable yipee
    dudz = Gu
    dTdz = GT
    drhodz = Grho
    dpdz = Gp
    dMdz = GM*M

    dYidz = (1 - M**2) * dYidx

    #residual check
    #res = Grho+Gu+(1-M**2)*(1/A)*dAdx
    #print(res)
    #print(str(rho*A*u))
    
    #h = DgasTPY.enthalpy_mass
    #print(h+(0.5*u**2))
    
    #print(str(u/(m.sqrt(gamma*R*T))))
    
    

    # output the subsonic state vector
    out6 = np.array([dxdz, dudz, dTdz, drhodz, dpdz, dMdz])
    #print(len(np.hstack((out6, dYidz))))
    return np.hstack((out6, dYidz))


# Subsonic function
def funSuperSonic(z, y):
    #unpack the state vector
    x = y[0]
    u = np.exp(y[1])
    T = np.exp(y[2])
    rho = np.exp(y[3])
    p = np.exp(y[4])
    M = y[5]
    YN = y[6:]  # These are mass fraction
    #print(M)


    # Decault gas for geeding to funcitons when one is needed
    DgasTPY = ct.Solution(IV.yaml)

    # We will need the vector later multiple times, so it gets computed here once
    #vector = makeVector(YN, DgasTPY)
    vector = YN

    # we need to calculate some standard gas properties from our system
    Wi = np.array(DgasTPY.molecular_weights)


    # mean molecular weight aswell
    #print([T, p, vector])
    DgasTPY.TPY = T, p, vector
    W = DgasTPY.mean_molecular_weight
    #print(W)


    # for reaction rates we will need a gas at the current temperature and pressure
    DgasTPY.TPY = T, p, vector
    widot = np.array(DgasTPY.net_production_rates)


    # current specific gas contant
    R = ct.gas_constant / W

    # specific heats at constant pressure, at temp and pressure, and ratio of specific heats
    cpi = np.array(DgasTPY.partial_molar_cp) / Wi
    cp = DgasTPY.cp_mass
    cv = DgasTPY.cv_mass
    gamma = cp / cv

    # enthalpies too
    hi = np.array(DgasTPY.partial_molar_enthalpies) / Wi

    # other length dependent funcitons we should only call once
    A = Af(x)
    qw = qwf(x)
    Pw = Pwf(x)
    dAdx = Derivative(x, Af, 10**(-3),)
    Cf = Cff(x)
    Dh = Dhf(x)

    # we need mdot aswell
    mdot = rho * A * u
    
    # The first thing that we need to do is compute the gas properties from species data
    # To start this we need the rate of change for every species per distance
    dYidx = dYidxf(u, rho, Wi, widot)

    # rest of the mixture props
    dWdx = dWdxf(YN, dYidx, Wi, W)
    dRdx = dRdxf(dWdx, W, R)
    
    
    
    #there is deviation from the primary source here, Cha states that this is a caloricaly perfect gas, however cantera no do that cp is temperature dependent
    #This means that we need to redetermine dcpdx assuming that cp is a function of temperature so we must do d/dx ourselves
    
    #because this is a reacting flow, enthalpy is not constant, so we need to compute the change in enthalpy over distance
    dh0dx = dh0dxf(mdot, qw, Pw)
    
    #this is also dependent on ho and sumterm so those go ehre too
    sumterm = sumtermf(hi, dh0dx, dYidx)

    #we need dTdx for this formulation
    GT = GTf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx)
    dTdx = T*(GT/(1-M**2))
    
    #we also need dcp/dT this needs to be found numerically
    dcpdT = Derivative(T, cpTf, 10**(-3),  p, vector, DgasTPY)
    dcpdx = dcpdxf(dYidx, cpi, dcpdT, dTdx)
    

    # ratio of specific heats type shit
    dgammadx = dgammadxf(T, p, gamma, cp, R, dcpdx, dRdx)

    # we now need to determine the G functions
    Gu = Guf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx)
    #GT calculated earlier
    Grho = Grhof(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx)
    Gp = Gpf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx)
    GM = GMf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx, dgammadx)

    # [x, u, T, rho, p, M, YN]
    # now we need to compute the derivatives of all our state variables with respect to z
    dxdz = abs(1 - M**2)
    
    #These variables were transformed into log form so that the solver is more stable yipee
    dudz = -Gu
    dTdz = -GT
    drhodz = -Grho
    dpdz = -Gp
    dMdz = -GM*M

    dYidz = (M**2 - 1) * dYidx

    #residual check
    #res = Grho+Gu+(1-M**2)*(1/A)*dAdx
    #print(res)
    #print(str(rho*A*u))
    
    #h = DgasTPY.enthalpy_mass
    #print(h+(0.5*u**2))
    
    #print(str(u/(m.sqrt(gamma*R*T))))
    
    

    # output the subsonic state vector
    out6 = np.array([dxdz, dudz, dTdz, drhodz, dpdz, dMdz])
    #print(len(np.hstack((out6, dYidz))))
    return np.hstack((out6, dYidz))


# This function produces a mass fraction vector for cantera
def makeVector(YN, Dgas):
    # create the result string
    result = ""

    # Go through all of the mass fractions and find what gas they are, then add them to the gas vector
    for pos, massfrac in enumerate(YN):
        result += Dgas.species_names[pos] + ":" + str(massfrac) + ", "

    # remove the extra comma and space
    result = result[:-2]

    return result


# function that gets the list of all species change rates over distance
def dYidxf(u, rho, Wi, widot):
    # equation 3 with mass addition term removed
    dYidx = (Wi * widot) / (rho * u)

    return dYidx


# function for determining the alteration in the mean molecular weight per distance
def dWdxf(YN, dYiDx, Wi, W):
    # compute mean molecular weight alteration over distance, equation 6
    dWdx = -((W) ** 2) * sum((1 / Wi) * dYiDx)

    return dWdx


# function for determining the chance of the specific gas constant per unit distance
def dRdxf(dWdx, W, R):
    # compute (1/R)dR/dx from equation 7
    dRdx = -(1 / W) * dWdx * R

    return dRdx

#function for determingin dcp/dt
def cpTf(T, p, vector, Dgas):
    #Tgets the value of the cp for some aribrary temperature
    Dgas.TPY = T, p, vector
    cp = Dgas.cp_mass
    
    return cp
 
# function for determining the change in specific heat of the total gas per unit distance
def dcpdxf(dYidx, cpi, dcpdT, dTdx):
    # calculate dxpdx from equation 8
    dcpdx = sum(cpi * dYidx)+(dcpdT*dTdx)

    return dcpdx


# function for determining the rate of change of gamma along distance
def dgammadxf(T, p, gamma, cp, R, dcpdx, dRdx):
    # calculate dgammadx from equation 9
    term1 = (-R / (cp * (cp - R))) * dcpdx
    term2 = (1 / (cp - R)) * dRdx
    dgammadx = gamma * (term1 + term2)

    return dgammadx


def dh0dxf(mdot, qw, Pw):
    # calculate dh0dx from equation equaiton 4
    dh0dx = (1 / mdot) * (-qw * Pw)
    #print(dh0dx)

    return dh0dx


# this is where we compute all of the "functions including forcing terms in the differential equation for flow peoperty phi"
def sumtermf(hi, dh0dx, dYidx):
    # this is the summation term frequently used in the other equations
    sumterm = dh0dx - sum(hi * dYidx)

    return sumterm


# [x, u, T, rho, p, M, YN]


# These are equations 17-21
def Guf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx):
    term1 = 1 / (cp * T)
    termR = (1 / R) * dRdx
    termA = -(1 / A) * dAdx
    termF1 = gamma * M**2
    termF2 = (2 * Cf) / Dh

    Gu = (term1 * sumterm) + termR + termA + (termF1 * termF2)
    return Gu


def GTf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx):
    term1 = (1 - (gamma * M**2)) / (cp * T)
    term2 = (gamma - 1) * M**2
    termR = -(1 / R) * dRdx
    termA = (1 / A) * dAdx
    termF1 = -gamma * M**2
    termF2 = (2 * Cf) / Dh

    GT = (term1 * sumterm) + (term2 * (termR + termA + (termF1 * termF2)))
    return GT


def Grhof(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx):
    term1 = -1 / (cp * T)
    termR = -(1 / R) * dRdx
    termM = M**2
    termA = (1 / A) * dAdx
    termF = -gamma * ((2 * Cf) / Dh)

    Grho = term1 * sumterm + termR + (termM * (termA + termF))
    return Grho


def Gpf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx):
    term1 = gamma * M**2
    term2 = (-1) / (cp * T)
    termR = (-1 / R) * dRdx
    termA = (1 / A) * dAdx
    term3 = -(1 + (gamma - 1) * M**2)
    termF = (2 * Cf) / Dh

    Gp = term1 * ((term2 * sumterm) + termR + termA + (term3 * termF))
    return Gp


def GMf(sumterm, T, M, cp, R, gamma, Cf, Dh, A, dRdx, dAdx, dgammadx):
    term1 = -(1 - M**2) / (2 * gamma)
    term2 = (1 + (gamma * M**2)) / 2
    term3 = 1 / (cp * T)
    termR = (1 / R) * dRdx
    term4 = -(1 + (((gamma - 1) / 2) * M**2))
    termA = (1 / A) * dAdx
    termF1 = -gamma * M**2
    termF2 = (2 * Cf) / Dh

    GM = (
        (term1 * dgammadx)
        + (term2 * ((term3 * sumterm) + termR))
        + (term4 * (termA + (termF1 * termF2)))
    )
    return GM


# Extranius fucntions that will need to be done differently
def qwf(x):
    return 0


def Pwf(fx):
    return 0


def Af(x):
    A = (m.pi * (RatL(x / 0.0254)) ** 2) * 0.0254**2
    #print(A)
    return A


def Cff(x):
    return 0


def Dhf(x):
    return 1


#sonic condition for choked flow
def sonic(z, y):
    return (y[5] - 1)
sonic.terminal = True


#length max for choked flow
def lengthThroat(z, y):
    return (y[0] - (LT*0.05+Lt)*0.0254)
lengthThroat.terminal = True


def lengthTotal(z, y):
    return (y[0] - LT*0.0254)
lengthTotal.terminal = True


#Function for shooting at a chamber pressure that prouces some specific mass flow #need to find the minimum Pc value for choked flow
def chamberpressureTMshooting(Pc, Tc, mdot, OF, yaml):
    #we need to determine density, velocity, and mach number
    g = ct.Solution(yaml)
    g = equalibrateHPWrapper(Pc, OF, 10**(-6))
    YN = g.Y
    g.TPY = Tc, Pc, YN
    
    #gas props
    gamma = g.cp_mass / g.cv_mass
    R = ct.gas_constant / g.mean_molecular_weight
    rho = g.density
    
    #we are going to need to produce a y0 vector
    x = 0
    u = mdot/(rho * Af(x))
    M = u/m.sqrt(gamma*R*Tc)
    
    y0 = np.array([x, np.log(u), np.log(Tc), np.log(rho), np.log(Pc), M])
    y0 = np.hstack((y0, np.array(YN)))
    
    #now that we have the y0 vector we must run the solver
    sol = solve_ivp(funSubsonic, [0, 10], y0, method="BDF", events=[sonic, lengthThroat], rtol=10**(-5), atol=10**(-6))
    
    #returns the solution currently
    return sol.y[5]


#functional version of chamberpressureTMshooting that takes in a vector for pc
def chamberpressureTMshootingParallel(Pcs, Tc, mdot, OF, yaml):
    #tun the function in paralell a bunch
    out = Parallel(n_jobs=IV.CC, backend="loky")(
        delayed(chamberpressureTMshooting)(Pc, Tc, mdot, OF, yaml) 
        for Pc in Pcs)
    
    #return all of the solutions to those differential equations
    return out


#uhh
def chamberpressureTMshootingParallelBisect(minPc, maxPc, Tc, mdot, OF, yaml, tol):
    #set error and tolerance, start the iterative loop
    error = 1
    while error > tol:
        #we need out list of pressures to compute at
        Pcs = np.linspace(minPc, maxPc, IV.CC)
        #print([minPc, maxPc])
        
        #we compute at these
        machs = chamberpressureTMshootingParallel(Pcs, Tc, mdot, OF, yaml)
        
        #we now must determine which values in there are the new bounds
        i = 0
        Found = False
        while (i < IV.CC-1) and not Found:
            #unpack the two solutions in question and the last mach number in them
            M1 = machs[i][-1]
            M2 = machs[i+1][-1]
            
            truth1 = M1 > (1-(10**(-5))) #choked condition
            truth2 = M2 < (1-(10**(-5))) #not choked
            
            #both of these must be true
            Found = truth1 and truth2
            #print([M1, M2])
            i+=1
        
        if not Found:
            print("Fuck")
        
        #now i is the index number for our machs, so we set our ounds to this new value
        minPc = Pcs[i-1]
        maxPc = Pcs[i]
        
        
        
        error = 1 - minPc/maxPc
        
    return([minPc, Pcs[i-2]])


def runCEA2():
    ceaOut = runCEA()
    #Starting conditions for iteration
    Pc: float = ceaOut[0].P
    Tc: float = ceaOut[0].T
    Tol: float = 10 ** (-3)
    RelError: float = 1
    Damp: float = 0.7
    

    Mdot: float = 0
    #Creates CombustionGas
    CombustionGas = ct.Solution(IV.yaml)
    #print(CombustionGas.report())

    i = 0
    while (RelError > Tol) and (i < 250):
        i += 1  # This is here to ensure no infinite loops
        PcOld: float = Pc #Keeps track of old Pc to find reletive error
        TcOld: float = Tc #Keeps track of old Tc to find reletive error


        # Gets Mdot for Both fuel and Ox sides all orifices, also total Mdot,
        FuelMdot: float = IV.FuelOrificeNum * Inj.MdotSPIONLY( IV.FuelOrificeCd, IV.FuelOrificeDiameter, IV.Fuel, IV.FuelTankT, Pc, IV.FuelTankP)
        BLCMdot: float = IV.BLCOrificeNum * Inj.MdotSPIONLY( IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, Pc, IV.FuelTankP)
        OxMdot: float = IV.OxOrificeNum * Inj.MdotSPIONLY(IV.OxOrificeCd, IV.OxOrificeDiameter, IV.Ox, IV.OxTankT, Pc, IV.OxTankP)
        Mdot = FuelMdot+OxMdot+BLCMdot
        Mdotcomb = FuelMdot+OxMdot
        #print(FuelMdot)
        #print(OxMdot)
        
        #Calculates OF ratio, technically not efficient to have it here or caculated this way, but eh.
        OF: float = OxMdot / FuelMdot

        #get the combustion gases at the current chamber pressure and OF ratio
        CombustionGas = equalibrateHPWrapper(Pc, OF, 10**(-6))
        
        #Calculates new ChamberPressure from the new tempature info, and gas propertys
        γ = CombustionGas.cp/CombustionGas.cv
        R = ct.gas_constant/CombustionGas.mean_molecular_weight
        Tc = CombustionGas.T
        #print([110000, Pc*0.5, Tc, Mdotcomb, OF, IV.yaml, 10**(-6)])
        Pc2 = chamberpressureTMshootingParallelBisect(Pc*0.8, Pc*1.2, Tc, Mdotcomb, OF, IV.yaml, 10**(-4))
        Pc = Pc2[0]
        #print((γ*R*Tc)**0.5)
        #print(Pc)
 
        #Calculates Max reletive error between ChamberPressure and ChamberTempature
        RelError = max((abs(Pc-PcOld))/(Pc),(abs(Tc-TcOld))/(Tc))
        #print("err:  " + str(RelError))

        #Calculates the chamge in ChamberPressure and ChamberTempature, makes sure its not so large it just overshoots everything and the engine "explodes"
        Pc -= Damp*(Pc-PcOld)
        Tc -= Damp*(Tc-TcOld)
 
 
    #we now need to use the converged solution to calculate the states everywhere
    #one final state vector :D
    #we need to determine density, velocity, and mach number
    g = ct.Solution(IV.yaml)
    g = equalibrateHPWrapper(Pc, OF, 10**(-6))
    YN = g.Y
    g.TPY = Tc, Pc, YN
    
    #gas props
    gamma = g.cp_mass / g.cv_mass
    R = ct.gas_constant / g.mean_molecular_weight
    rho = g.density
    
    #we are going to need to produce a y0 vector
    x = 0
    u = Mdot/(rho * Af(x))
    M = u/m.sqrt(gamma*R*Tc)
    
    y0 = np.array([x, np.log(u), np.log(Tc), np.log(rho), np.log(Pc), M])
    y0 = np.hstack((y0, np.array(YN)))
    
    #now that we have the y0 vector we must run the solver
    solsub = solve_ivp(funSubsonic, [0, 10], y0, method="BDF", events=[sonic, lengthThroat], rtol=10**(-5), atol=10**(-6))
    
    
    
    #I LIED ONE MORE STATE VECTOR AAAAAAAAAAAAA ITS 2 MORE
    #we need the previous unsolved choked case.
    x = 0
    u = Mdot/(rho * Af(x))
    M = u/m.sqrt(gamma*R*Tc)
    
    y0 = np.array([x, np.log(u), np.log(Tc), np.log(rho), np.log(Pc2[1]), M])
    y0 = np.hstack((y0, np.array(YN)))
    
    #now that we have the y0 vector we must run the solver
    solsub2 = solve_ivp(funSubsonic, [0, 10], y0, method="BDF", events=[sonic, lengthThroat], rtol=10**(-5), atol=10**(-6))
    
    dety0R = dety0(solsub, solsub2)
    y02 = dety0R[0]
    #print(solsub.y[0, -1])
    #print(y02[0])
    #y02[0] = solsub.y[0, -1]
    #we gonna make mach 1
    #print(y02)
    
    
    solsup = solve_ivp(funSuperSonic, [0, 10], y02, method="BDF", events=lengthTotal, rtol=10**(-6), atol=10**(-7))
    
    
    #as this is our final subsonic state we need to under the log on u temp pressure and density
    solsub.y[1:5] = np.exp(solsub.y[1:5])
    solsup.y[1:5] = np.exp(solsup.y[1:5])
    return CombustionGas, Mdot, solsub, solsup, dety0R[1]


#function for determining the initial state of y0 for the supersonic solution
def dety0(sol, solprev):
    #we need the derivative of the last state vector and some data from it
    #these state vectors are loged
    yk = sol.y[:, -1].copy()  #get last state vector
    ykm1 = solprev.y[:, -1].copy() #the one before that
    
    #compute the derivative of these
    dykdz = funSubsonic(0, yk)
    dykm1dz = funSubsonic(0, ykm1)
    #these have dMdz = GM*M
    
    #we need both GMS
    GMk = dykdz[5]/yk[5]
    GMkm1 = dykm1dz[5]/ykm1[5]
    
    #this formulation will also need unlogged ys and now that we ve computed the derivatives this is fine to do
    yk[1:5] = np.exp(yk[1:5])
    ykm1[1:5] = np.exp(ykm1[1:5])
    
    #equation45
    ystar = yk + ((yk-ykm1)/(GMk-GMkm1))*(0-GMk)
    
    #we need several more bits, Mm1 M and dydM
    #dydM = (dy/dz)/(dM/dz)
    Mc = 1+10**(-4)
    
    #mk here needs to be the point when M is 0.999
    index = np.argmin(np.abs(sol.y[5] - 0.999))
    yminus = sol.y[:, index].copy()
    Mk = yminus[5]
    
    dMdz = GMk*Mk
    dydM = dykdz/dMdz
    
    #now we need to compute the a and b functions eq 49 and 50
    a = (1/(1-Mk))*(((ystar-yminus)/(1-Mk))-dydM)
    b = dydM - 2*((ystar-yminus)/(1-Mk))
    
    #equation 48
    def yIf(M): 
        return a*(1-M)**2 + b*(1-M) + ystar
    y0 = yIf(Mc)
    #print(y0)
    #we need to log this again
    y0[1:5] = np.log(y0[1:5])
    
    #while we are here we also need the data to fill this interpolated point
    print((y0[0]-yk[0])/0.0254)
    #print(y0)
    print(Mk)
    machs = np.linspace(Mk, Mc, 100)
    states = np.stack([yIf(M) for M in machs], axis=1)
    
    return [y0, states]



