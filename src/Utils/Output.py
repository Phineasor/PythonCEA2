import cea as cea
import cantera as ct

ceaOut = cea.runCEA()

γ = ceaOut[0].cp/ceaOut[0].cv
R = ct.gas_constant/ceaOut[0].mean_molecular_weight

print(ceaOut[0]())
#print(R)