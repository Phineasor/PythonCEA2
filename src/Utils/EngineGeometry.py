# Calculates the needed engine geometry paramaters based on the given infromation in the input file
"""
@author: phineas
"""

import math
from ..Inputs.InputValues import(
    ExpRatio,
    ConRatio,
    Dt,
    Rcont,
    ThetaCont,
    Rexp,
    ThetaExp,
    Lchamber,
    Method,
)

d2r = math.pi/180

if Method == 1:
    # Calculates needed Engine dimension values from Grissom
    Dc = Dt * math.sqrt(ConRatio)
    De = Dt * math.sqrt(ExpRatio)
    D2 = Dc - 2 * Rcont * (1 - math.cos(ThetaCont*d2r))
    D3 = Dt + 2 * Rexp * (1 - math.cos(ThetaCont*d2r))
    D5 = Dt + 2 * Rexp * (1 - math.cos(ThetaExp*d2r))
    L2 = Lchamber + Rcont * math.sin(ThetaCont*d2r)
    L3 = L2 + ((D2 - D3) / (2 * math.tan(ThetaCont)))
    Lt = L3 + Rexp * math.sin(ThetaCont*d2r)
    L5 = Lt + Rexp * math.sin(ThetaExp*d2r)
    LT = L5 + ((De - D5) / (2 * math.tan(ThetaExp*d2r)))
    # End of engine dimension value calculation


def RatL(L):
    # Start of contour calculation
    if L < Lchamber:
        return 0.5 * (Dc)
    if L < L2:
        return 0.5 * (Dc - 2 * (Rcont - math.sqrt(Rcont**2 - (L - Lchamber) ** 2)))
    if L < L3:
        return 0.5 * (D2 - 2 * (L - L2) * (math.tan(ThetaCont*d2r)))
    if L < L5:
        return 0.5 * (Dt + 2 * (Rexp - math.sqrt(Rexp**2 - (L - Lt) ** 2)))
    else:
        return 0.5 * (D5 + 2 * (L - L5) * (math.tan(ThetaExp*d2r)))
