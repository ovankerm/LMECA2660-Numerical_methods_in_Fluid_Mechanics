import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl

for i in range(1, 11):
    t = 1 * (11-i)

    fig, ax = plt.subplots(1, 3, figsize=(14, 7))

    fig.suptitle("t = %d"%t)

    with open("Project/output/T-%d.txt"%t, 'r') as f:
        l = f.readline().split()
        Nx = int(l[2])
        Ny = int(l[5])
        h = float(l[8])
        T = np.zeros((Ny-1, Nx-1))
        for i in range((Nx-1) * (Ny-1)):
            l = f.readline().split()
            T[int(l[5]) - 1, int(l[2]) - 1] = float(l[8])

    X_T, Y_T = np.meshgrid(np.linspace(h, 2/3 - h, endpoint=True, num=(Nx-1)), np.linspace(h, 1 - h, endpoint=True, num=(Ny-1)))

    cs_T = ax[0].contourf(X_T, Y_T, T, cmap='bwr', levels=100)
    ax[0].set_title("T")
    plt.colorbar(cs_T, ax=ax[0])

    with open("Project/output/V-%d.txt"%t, 'r') as f:
        l = f.readline().split()
        Nx = int(l[2])
        Ny = int(l[5])
        h = float(l[8])
        U = np.zeros((Ny, Nx-1))
        for i in range((Nx-1) * (Ny)):
            l = f.readline().split()
            U[int(l[5]), int(l[2]) - 1] = float(l[8])

    X_U, Y_U = np.meshgrid(np.linspace(h, 2/3-h, endpoint=True, num=(Nx-1)), np.linspace(0, 1, endpoint=True, num=(Ny)))

    cs_U = ax[1].contourf(X_U, Y_U, U, cmap='bwr', levels=100)
    ax[1].set_title("V")
    plt.colorbar(cs_U, ax=ax[1])

    # with open("Project/output/0_V-%d.txt"%t, 'r') as f:
    #     l = f.readline().split()
    #     Nx = int(l[2])
    #     Ny = int(l[5])
    #     h = float(l[8])
    #     U = np.zeros((Ny, Nx-1))
    #     for i in range((Nx-1) * (Ny)):
    #         l = f.readline().split()
    #         U[int(l[5]), int(l[2]) - 1] = float(l[8])

    # X_U, Y_U = np.meshgrid(np.linspace(h, 2/3-h, endpoint=True, num=(Nx-1)), np.linspace(0, 1, endpoint=True, num=(Ny)))

    # cs_U = ax[1].contourf(X_U, Y_U, U, cmap='bwr', levels=100)
    # ax[1].set_title("0_V")
    # plt.colorbar(cs_U, ax=ax[1])

    with open("Project/output/U-%d.txt"%t, 'r') as f:
        l = f.readline().split()
        Nx = int(l[2])
        Ny = int(l[5])
        h = float(l[8])
        V = np.zeros((Ny, Nx-1))
        for i in range((Nx-1) * (Ny)):
            l = f.readline().split()
            V[int(l[5]), int(l[2]) - 1] = float(l[8])

    X_V, Y_V = np.meshgrid(np.linspace(h, 2/3 - h, endpoint=True, num=(Nx-1)), np.linspace(0, 1, endpoint=True, num=Ny))

    cs_V = ax[2].contourf(X_V, Y_V, V, cmap='bwr', levels=100)
    ax[2].set_title("U")
    plt.colorbar(cs_V, ax=ax[2])


    # fig, ax = plt.subplots(1, 3, figsize=(14, 7))

    # fig.suptitle("t = %.2f"%t)

    # with open("Project/output/P-%.4f.txt"%t, 'r') as f:
    #     l = f.readline().split()
    #     Nx = int(l[2])
    #     Ny = int(l[5])
    #     h = float(l[8])
    #     T = np.zeros((Ny-1, Nx-1))
    #     for i in range((Nx-1) * (Ny-1)):
    #         l = f.readline().split()
    #         T[int(l[5]), int(l[2])] = float(l[8])

    # X_T, Y_T = np.meshgrid(np.linspace(h, 2/3 - h, endpoint=True, num=(Nx-1)), np.linspace(h, 1 - h, endpoint=True, num=(Ny-1)))

    # cs_T = ax[0].contourf(X_T, Y_T, T, cmap='bwr', levels=100)
    # ax[0].set_title("P")
    # plt.colorbar(cs_T, ax=ax[0])

    # with open("Project/output/U_star-%.4f.txt"%t, 'r') as f:
    #     l = f.readline().split()
    #     Nx = int(l[2])
    #     Ny = int(l[5])
    #     h = float(l[8])
    #     U = np.zeros((Ny-1, Nx))
    #     for i in range((Nx) * (Ny-1)):
    #         l = f.readline().split()
    #         U[int(l[5]) - 1, int(l[2])] = float(l[8])

    # X_U, Y_U = np.meshgrid(np.linspace(0, 2/3, endpoint=True, num=Nx), np.linspace(h, 1 - h, endpoint=True, num=(Ny-1)))

    # cs_U = ax[1].contourf(X_U, Y_U, U, cmap='bwr', levels=100)
    # ax[1].set_title("U")
    # plt.colorbar(cs_U, ax=ax[1])

    # with open("Project/output/V_star-%.4f.txt"%t, 'r') as f:
    #     l = f.readline().split()
    #     Nx = int(l[2])
    #     Ny = int(l[5])
    #     h = float(l[8])
    #     V = np.zeros((Ny, Nx-1))
    #     for i in range((Nx-1) * (Ny)):
    #         l = f.readline().split()
    #         V[int(l[5]), int(l[2]) - 1] = float(l[8])

    # X_V, Y_V = np.meshgrid(np.linspace(h, 2/3 - h, endpoint=True, num=(Nx-1)), np.linspace(0, 1, endpoint=True, num=Ny))

    # cs_V = ax[2].contourf(X_V, Y_V, V, cmap='bwr', levels=100)
    # ax[2].set_title("V")
    # plt.colorbar(cs_V, ax=ax[2])

plt.show()