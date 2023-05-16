import numpy as np
from matplotlib import pyplot as plt

N = 1000

Tavg = np.zeros(N)
sigT = np.zeros(N)
iter = np.zeros(N)
qe = np.zeros(N)
Ts = np.zeros(N)
Rehw = np.zeros(N)
Reh = np.zeros(N)

with open("Project/output/diagnostics.txt", "r") as f:
    for i in range(N):
        l = f.readline().split()
        Tavg[i] = float(l[5])
        sigT[i] = float(l[8])
        iter[i] = float(l[2])
        qe[i] = float(l[11])
        Ts[i] = float(l[14])
        Rehw[i] = float(l[17])
        Reh[i] = float(l[20])

plt.figure()
plt.plot(iter, sigT, label='sigT')
plt.legend()
plt.figure()
plt.plot(iter, Tavg, label='Tavg')
plt.legend()
plt.figure()
plt.plot(iter, qe, label='qe')
plt.legend()
plt.figure()
plt.plot(iter, Ts, label='Ts')
plt.legend()
plt.figure()
plt.plot(iter, Rehw, label='Rehw')
plt.legend()
plt.figure()
plt.plot(iter, Reh, label='Reh')
plt.legend()
plt.show()