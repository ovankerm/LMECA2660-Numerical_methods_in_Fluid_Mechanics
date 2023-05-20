import numpy as np
from matplotlib import pyplot as plt

iter = 1572000
N = 15
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

y = np.linspace(0, 1, num=Ny, endpoint=True)

ind = int(0.7 * Ny)
x = np.linspace(-h/2, -h/2 + (N-1) * h, num = N, endpoint=True)

fig, ax = plt.subplots(figsize=(9,6))
ax.grid()
for i in range(1):
    ax.plot(v[ind + 7 * i, :], x, 'k', label=f'y/H = {y[ind + 7 * i]}')
ax.set_ybound(0)
ax.legend()
ax.set_xlabel(r'$\frac{v}{U}$', fontsize=15)
ax.set_ylabel(r'$\frac{x}{H}$', fontsize=15)
ax.set_title('Vertical component of the velocity near the wall')

fig.savefig('Project/images/boundary_layer.eps', format='eps')
plt.show()

