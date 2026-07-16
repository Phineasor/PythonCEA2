#This class provides functions that convert lists from propellent names into
#A) Coolprop valid chemicals
#B) Returns some spefic cantera fuel/Y vector

#input tables and output tables
Input1 = ["Ethanol", "Kerosene", "LOx"] #Oxidizers need to be at the end here
Input2 = ["LOx"]
OutputCool = ["Ethanol", "nDodecane", "O2"]

#Function A
def conCool(propellent):    
    return OutputCool[Input1.index(propellent)]
#https://chemistry.cerfacs.fr/en/chemical-database/mechanisms-list/dagauts-mechanism/
#Function B
def conCant(fuel, ox, OF):
    OutputCant = [["O2:"+str(OF)+", C2H5OH:1"],                                #Ethanol mixture with: LOx, .....
                  ["O2:"+str(OF)+", NC10H22:0.74, PHC3H7:0.15, CYC9H18:0.11"]] #Kerosene mixture with: LOx, .....
    return OutputCant[Input1.index(fuel)][Input2.index(ox)]