import numpy as np
from matplotlib import pyplot as plt

N = 400000

Tavg_no = np.zeros(N)
sigT_no = np.zeros(N)
iter_no = np.zeros(N)

with open("Project/output/diagnostics_no.txt", "r") as f:
    for i in range(N):
        l = f.readline().split()
        Tavg_no[i] = float(l[5])
        sigT_no[i] = float(l[8])
        iter_no[i] = i

Tavg = np.zeros(N)
sigT = np.zeros(N)
iter = np.zeros(N)

with open("Project/output/diagnostics.txt", "r") as f:
    for i in range(N):
        l = f.readline().split()
        Tavg[i] = float(l[5])
        sigT[i] = float(l[8])
        iter[i] = i

plt.figure()
plt.plot(iter_no, sigT_no)
plt.plot(iter, sigT)
plt.show()