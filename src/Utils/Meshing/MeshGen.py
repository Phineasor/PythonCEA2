import math
import numpy as np
import matplotlib.pyplot as plt
from ..Bisect import Bisect
from ...Inputs import IVMesh as M
from ...Inputs import InputValues as IV


points = np.zeros((IV.CellNum, M.BLnum+M.FFres))

x = np.append(np.linspace(3, 3.125, 20), np.delete(np.linspace(3.125, 5, 31), 0))
y = [0]*50
print(points)


for i in range(points):
    for j in range(points[0]):
        points[i][j] = 


plt.rcParams['figure.dpi'] = 500
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111)
ax.set_facecolor('black')

ax.scatter(x, y, color='white',linestyle='--', linewidth=0.1)

plt.title("3D engine model")
plt.legend(fancybox=False, shadow=True, framealpha=1,fontsize='small',loc='lower left')
plt.savefig("test2")
