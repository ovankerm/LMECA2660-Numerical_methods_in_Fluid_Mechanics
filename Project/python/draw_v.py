import numpy as np
from matplotlib import pyplot as plt
fig, ax = plt.subplots(1, 2, figsize=(13, 6))

with open("Project/output/U_star.txt", 'r') as f:
    Nx = 8
    Ny = 8
    h = 0.1
    V = np.zeros((Ny-1, Nx))
    for i in range((Ny-1) * (Nx)):
        l = f.readline().split()
        V[int(l[5]) - 1, int(l[2])] = float(l[8])
    X_V, Y_V = np.meshgrid(np.linspace(0, (Nx-1) * h, endpoint=True, num=Nx), np.linspace(h/2, (Ny-1) * h - h/2, endpoint=True, num=(Ny-1)))
    cs_V = ax[0].contourf(X_V, Y_V, V, cmap='bwr', levels=100)
    ax[0].set_title("U")
    plt.colorbar(cs_V, ax=ax[0])
    ax[0].set_xlim(0, 0.7)
    ax[0].set_ylim(0, 0.7)

with open("Project/output/V_star.txt", 'r') as f:
    Nx = 8
    Ny = 8
    h = 0.1
    V = np.zeros((Ny, Nx-1))
    for i in range((Nx-1) * (Ny)):
        l = f.readline().split()
        V[int(l[5]), int(l[2]) - 1] = float(l[8])
    X_V, Y_V = np.meshgrid(np.linspace(h/2, (Nx-1) * h - h/2, endpoint=True, num=(Nx-1)), np.linspace(0, (Ny-1) * h, endpoint=True, num=(Ny)))
    cs_V = ax[1].contourf(X_V, Y_V, V, cmap='bwr', levels=100)
    ax[1].set_title("V")
    plt.colorbar(cs_V, ax=ax[1])
    ax[1].set_xlim(0, 0.7)
    ax[1].set_ylim(0, 0.7)

plt.show()