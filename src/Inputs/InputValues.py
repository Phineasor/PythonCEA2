# General InputValues
yaml = "PCRL-Mech1.yaml"


# Input numbers in the following section are in inches and degrees(expet the ratios they are ratios)
ExpRatio = 3.291
ConRatio = 9
Dt = 1
Rcont = 0.5
ThetaCont = 45
Rexp = 0.5
ThetaExp = 15
Lchamber = 3
# End of input section for chamber and nozzle geometry

# Method of radius calculation (linear/parabolic, 1/2)
Method = 1
# End of Method input

# Number of cells in the 1D CEA
CellNum = 250
# End of cells

# Injector charataristics, Ox holes, fuel holes BLC holes
FuelOrificeDiameter = 0.0550
FuelOrificeNum = 8
FuelOrificeCd = 0.7
FuelAngle = 0  # NOT CURRENTLY USED

OxOrificeDiameter = 0.0591
OxOrificeNum = 8
OxOrificeCd = 0.7
OxAngle = 0  # NOT CURRENTLY USED

BLCOrificeDiameter = 0.021 #0.0197
BLCOrificeNum = 16 * 0.5
BLCOrificeCd = 0.7
BLCAngle = 45  # I was wrong this is used now
# End of Injector pro

# Tank charataristics
FuelTankP = 350  # PSI
FuelTankT = 290  # K
Fuel = "Ethanol"

OxTankP = 350  # PSI
OxTankT = 85  # K
Ox = "O2"
# End Of tank charataristics

# Ambiant Conditions
AmbP = 101325  # In pascales, sorry not PSI here
AmbT = 300  # Ambient temp in K

# Film Cooling paramaters
Aw = 0.5  # Wall Absotrptivity
et = 0.2  # Estamation of turbulence correction factor

#Ray comp values
mainName = "AbsCoefData" #The on disk name should be mainName + CellNum
Memmap = True #should mem map be used, or load into main array
Xlocation = 0 #location to compute x at.

#true AbsCoefName
AbsCoefName = mainName + str(CellNum) + '.npy'
