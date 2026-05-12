# File to generate radiative power values
# fmt: off
"""
@author: phineas
"""
#External modules
import numpy as np
import math as m
from multiprocessing import shared_memory

#Internal modules
from cea import *
import InputValues as IV
from RayMarch import getRay
from concurrent.futures import ProcessPoolExecutor


n_cores = 24  # how many processes you want

inputs = np.array([(0, 0, 0)])
M = np.linspace(-90, 90, 35)
N = np.linspace(0, 90, 18)
M[0] = -89.99
M[-1] = 89.99
N[-1] = 89.99
#print(M)

i = 0
while i < len(M):
    j = 0
    while j < len(N):
        inputs = np.append(inputs, np.array([(IV.Xlocation, M[i]*(m.pi/180), N[j]*(m.pi/180))]), axis = 0)
        j+=1
    i+=1
inputs = np.delete(inputs, 0, 0)

def call_CompRay(args):
    return CompRay(*args)

#will use memory map to compute the rays
if (__name__ == "__main__") and IV.Memmap:
    from ComputeRays import CompRay
    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        results = list(executor.map(call_CompRay, inputs))
    results = np.array(results, dtype=np.float64)
    name = str(IV.Xlocation) + "PowerSter"
    np.save(name, results)

#will use shared ram memory to ocmpute the rays
if (__name__ == "__main__") and not IV.Memmap:
    AbsCoefArray = np.load(IV.AbsCoefName, allow_pickle=True)
    shm = shared_memory.SharedMemory(name = 'AbsCoefDataMemoryBuffer', create=True, size=AbsCoefArray.nbytes)
    from ComputeRays import CompRay

    with ProcessPoolExecutor(max_workers=n_cores) as executor:
        results = list(executor.map(call_CompRay, inputs))
    results = np.array(results, dtype=np.float64)
    name = str(IV.Xlocation) + "PowerSter"
    np.save(name, results)
    
if not (__name__ == "__main__") and not IV.Memmap:
    from ComputeRays import CompRay