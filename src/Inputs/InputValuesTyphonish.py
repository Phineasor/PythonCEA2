# General InputValues
yaml = "Dagaut_HP.yaml"


#Propellent charictarization graph generaion options
LongGraphs = False
ShortGraphs = False


#is this method allowed to run
Physical = True

# Input numbers in the following section are in inches and degrees(expet the ratios they are ratios)
ExpRatio = 5.85
ConRatio = 14.0625
Dt = 1.2
Rcont = 0.75
ThetaCont = 35
Rexp = 0.5
ThetaExp = 15
Lchamber = 5
# End of input section for chamber and nozzle geometry

# CEA Methods
Method = 1  #Method of radius calculation (linear/parabolic, 1/2)
cStarEff = 1
# End of Method input

# Number of cells in the 1D CEA
CellNum = 250
# End of cells

# Injector charataristics, Ox holes, fuel holes BLC holes
FuelOrificeDiameter = 0.0520
FuelOrificeNum = 16
FuelOrificeCd = 0.7
FuelAngle = 0  # NOT CURRENTLY USED

OxOrificeDiameter = 0.0670
OxOrificeNum = 16
OxOrificeCd = 0.7
OxAngle = 0  # NOT CURRENTLY USED

BLCOrificeDiameter = 0.021 #0.0197
BLCOrificeNum = 0
BLCOrificeCd = 0.7
BLCAngle = 45  # I was wrong this is used now
# End of Injector pro

# Tank charataristics
FuelTankP = 750  # PSI
FuelTankT = 290  # K
Fuel = "Kerosene"

OxTankP = 750  # PSI
OxTankT = 85  # K
Ox = "LOx"
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

#Normal Pc Thrust driven solver inputs to determin chamber values.
Dc = 4.5 #Arbitrary chamber diameter in inches
Pc = 4200000 #Chamber pressure in pascal
F = 5000 #engine thrust in newtons
OF = 2.18 #OF ratio
LStar = 60 #ideal Lstar in inches