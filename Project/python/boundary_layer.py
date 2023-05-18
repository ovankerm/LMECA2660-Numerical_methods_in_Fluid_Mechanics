import numpy as np
from matplotlib import pyplot as plt

iter = 1562568
N = 20
with open(f"Project/output/v-{iter}.txt", 'r') as f:
    l = f.readline().split()
    Nx = int(l[2])
    Ny = int(l[5])
    h = float(l[8])
    v = np.zeros((Ny, N))
    for i in range(Ny):
        for j in range(N):
            l = f.readline().split()
            v[int(l[1]), int(l[0])] = float(l[2])

ind = int(0.3 * Ny)
x = np.linspace(-h/2, -h/2 + (N-1) * h, num = N, endpoint=True)

fig, ax = plt.subplots()
ax.grid()
for i in range(10):
    ax.plot(v[ind + 2 * i, :], x)
ax.set_ybound(0)
plt.show()

