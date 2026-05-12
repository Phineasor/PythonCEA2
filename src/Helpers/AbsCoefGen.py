# File to generate absorbtion coefficient data from hapi
# fmt: off
"""
@author: phineas
"""
#External modules
import matplotlib.pyplot as plt
import numpy as np
from hapi import *
#from tables import *

#Internal modules
from cea import *
import InputValues as IV

#Tells Hapi where the data is stored
db_begin('Data')

#conversion factors
pa2atm = 9.86923*(10**(-6))

#Needed engine data to produce abscoeff
cea = runCEA()
AxVal = AxialValues(cea[0].T, cea[0].P, cea[0].density, cea[0])
#print(AbsCoefAray[240])

AbsCoefArayBad = np.array([None]*IV.CellNum)


i_H2O = cea[0].species_index('H2O')
i_CO  = cea[0].species_index('CO')
i_CO2 = cea[0].species_index('CO2')

moleFrac = cea[0].X[[i_H2O, i_CO, i_CO2]]
for i in range(IV.CellNum):
        nu,coef = absorptionCoefficient_HT(SourceTables=['H2O', 'CO', 'CO2'], WavenumberStep = 0.001, Environment = {'p':((AxVal[1][i])*pa2atm), 'T':AxVal[0][i]}, OmegaWingHW = 100, HITRAN_units = False, Diluent = {'CO2':moleFrac[0], 'CO':moleFrac[1], 'H2O':moleFrac[2]})
        AbsCoefArayBad[i] = [nu, coef]
        i+=1

AbsCoefAray3 = np.array([[[0.0]*len(AbsCoefArayBad[0][0]), [0.0]*len(AbsCoefArayBad[0][0])]]*IV.CellNum)
count = 0
i = 0
while i < IV.CellNum:
        j = 0
        while j < 2:
                k = 0
                while k < len(AbsCoefArayBad[0][0]):
                        AbsCoefAray3[i][j][k] = AbsCoefArayBad[i][j][k]
                        k += 1
                        count += 1
                        
                        if count % 100000 == 0:
                                print(count)
                j += 1
        i += 1
np.save('AbsCoefAray3', AbsCoefAray3)