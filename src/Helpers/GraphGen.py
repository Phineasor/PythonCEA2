#Produces many propellent charictarating graphs dependent on the input propellents, and selected theochemical properties.
#external imports
import math as m
import numpy as np
import cantera as ct

import CoolProp.CoolProp as CP
CP.set_config_bool(CP.ENABLE_SUPERANCILLARIES, True)

import matplotlib.pyplot as plt

#internal imports
from ..Inputs import InputValues as IV
from ..Helpers.PropName import conCool, conCant
from ..Utils.Bisect import Bisect, Derivative
from ..Utils.cea import equalibrateHPWrapper

#function for determining Tc and Isp at sea level from Pc and OF is dependent on true for Isp
#uses the prop qualities from the input values table
def getTempAndIsp(OF: float, Pc: float, IspOrTc: bool, tol: float):
    #gets the combustion gas results from the input values
    CombustionGas = equalibrateHPWrapper(Pc, OF, tol)
    
    #We must now compute the ideal specific impulse
    y = CombustionGas.cp/CombustionGas.cv
    R = ct.gas_constant/CombustionGas.mean_molecular_weight
    Tc = CombustionGas.T
    
    #Ideal characteristic velocity
    exp = (y+1)/(y-1)
    cStar = (m.sqrt(y*R*Tc))/(y*m.sqrt((2/(y+1))**exp))
    
    #Thrust Coeff: consider, the nozzle is iedally expanded, so  p2-p3=0 in sutton:3-30
    frac1 = (2*y**2)/(y-1)
    frac2 = (2/(y+1))
    frac3 = (y+1)/(y-1)
    CF = m.sqrt(frac1*(frac2**frac3)*(1-(101325/Pc)**((y-1)/y)))
    
    #specific impulse from cStar and CF
    Isp = (cStar*CF)/9.81
    
    #Changes output to isp if IspOrTc is true
    if IspOrTc:
        output = Isp
    else:
        output = Tc
    #print(output)
    return output



if IV.LongGraphs:
    #Production for Tensor of spefic impulse with respect to Pc and OF
    #Prodcution for Tensor of Isp with respect to Pc and OF
    Pc = np.linspace(500000, 5000000, 10)
    OF = np.linspace(0.5, 5, 50)

    VgetTempAndIsp = np.vectorize(getTempAndIsp)

    Ispfig = plt.figure()
    ax = Ispfig.add_subplot(projection='3d')

    PcMG, OFMG = np.meshgrid(Pc, OF)
    IspList = VgetTempAndIsp(OFMG, PcMG, True, 10**(-3))
    ax.plot_surface(OFMG, PcMG/1000, IspList, cmap="turbo")

    # Tweak the limits and add latex math labels.
    ax.set_xlabel("OF Ratio")
    ax.set_ylabel("Chamber Pressure [kPa]")
    ax.set_zlabel("Specific Impulse [s]")
    ax.view_init(elev=30, azim=-60 - 90)

    #Prodcution for Tensor of Tc with respect to Pc and OF
    Tempfig = plt.figure()
    ax = Tempfig.add_subplot(projection='3d')

    PcMG, OFMG = np.meshgrid(Pc, OF)
    TcList = VgetTempAndIsp(OFMG, PcMG, False, 10**(-3))
    ax.plot_surface(OFMG, PcMG/1000, TcList, cmap="turbo")

    #Tweak the limits and add latex math labels.
    ax.set_xlabel("OF Ratio")
    ax.set_ylabel("Chamber Pressure [kPa]")
    ax.set_zlabel("Combustion Temperature [K]")
    ax.view_init(elev=30, azim=-60 - 90)
    plt.show()



#A 2D verison of both of these plots at some set chamber pressure
if IV.ShortGraphs:
    OF = np.linspace(0.25, 5, 100)
    VgetTempAndIsp = np.vectorize(getTempAndIsp)

    IspList = VgetTempAndIsp(OF, 3500000, True, 10**(-3))
    TcList  = VgetTempAndIsp(OF, 3500000, False, 10**(-3))
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    
    # Tweak the limits and add latex math labels.
    ax1.set_title('Pressure: 3500 kPa')
    ax1.set_xlabel("OF Ratio")
    
    ax1.plot(OF, TcList,  color="orangered")
    ax1.set_ylabel("Chamber Temperature [K]", color="orangered")
    ax1.tick_params(axis="y", labelcolor="orangered")
    
    ax2.plot(OF, IspList, color="black")
    ax2.set_ylabel("Specific Impulse [s]", color="black")
    ax2.tick_params(axis="y", labelcolor="black")
    
    ax2.set_ylim(125, 300)
    ax1.set_ylim(850, 4000)
    ax2.set_yticks(np.linspace(125, 300, 11))
    ax1.set_yticks(np.linspace(850, 4000, 11))
    
    
    ax1.grid()
    ax2.grid()
    fig.tight_layout()
    plt.show()
 
 
def func(x, *args):
    return Derivative(getTempAndIsp, x, *args)
    
#print(Bisect(Derivative, 0.1, 5, 10**(-7), getTempAndIsp, 4200000, True, 10**(-8)))




