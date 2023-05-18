import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

levels_omega = np.linspace(-2, 2, num=50)
levels_T = np.linspace(0, 0.025, num=50)
levels_U = np.linspace(0, 0.05, num=50)
levels_Rehw = np.linspace(0, 40, num=50)
levels_Reh = np.linspace(0, 25, num=50)

size = (5,7.5)
diter = 1310
n_iter = 1000

fig, ax = plt.subplots(figsize=size)
theta = np.linspace(0, 2 * np.pi, num = 100)

def animate(t):
    ax.clear()

    with open("Project/output/T-%d.txt"%(diter*(t+1)), 'r') as f:
        l = f.readline().split()
        Nx = int(l[2])
        Ny = int(l[5])
        h = float(l[8])
        t = float(l[11])
        r = 1/5 * np.cos(3 * (theta - 0.1 * t))
        T = np.zeros((Ny-1, Nx-1))
        ax.set_title(f"t = {t}")
        for i in range((Nx-1) * (Ny-1)):
            l = f.readline().split()
            T[int(l[1]) - 1, int(l[0]) - 1] = float(l[2])


    X_T, Y_T = np.meshgrid(np.linspace(h/2, 2/3 - h/2, endpoint=True, num=(Nx-1)), np.linspace(h/2, 1 - h/2, endpoint=True, num=(Ny-1)))

    cs_T = ax.contourf(X_T, Y_T, T, cmap='inferno', levels=levels_T, extend='both')
    ax.fill(r * np.cos(theta) + 1/3, r * np.sin(theta) + 1/3, color='grey')
    ax.fill(1/25 * np.cos(theta) + 1/3, 1/25 * np.sin(theta) + 1/3, color='grey')

    # with open("Project/output/w-%d.txt"%(diter*(t+1)), 'r') as f:
    #     l = f.readline().split()
    #     Nx = int(l[2])
    #     Ny = int(l[5])
    #     h = float(l[8])
    #     t = float(l[11])
    #     r = 1/5 * np.cos(3 * (theta - 0.1 * t))
    #     T = np.zeros((Ny, Nx))
    #     ax.set_title(f"t = {t}")
    #     for i in range((Nx) * (Ny)):
    #         l = f.readline().split()
    #         T[int(l[1]), int(l[0])] = float(l[2])

    # X_T, Y_T = np.meshgrid(np.linspace(0, 2/3, endpoint=True, num=(Nx)), np.linspace(0, 1, endpoint=True, num=(Ny)))

    # cs_T = ax.contourf(X_T, Y_T, T, cmap='bwr', levels=levels_omega, extend='both')
    # ax.fill(r * np.cos(theta) + 1/3, r * np.sin(theta) + 1/3, color='grey')
    # ax.fill(1/25 * np.cos(theta) + 1/3, 1/25 * np.sin(theta) + 1/3, color='grey')

    # with open("Project/output/Rehw-%d.txt"%(diter*(t+1)), 'r') as f:
    #     l = f.readline().split()
    #     Nx = int(l[2])
    #     Ny = int(l[5])
    #     h = float(l[8])
    #     t = float(l[11])
    #     r = 1/5 * np.cos(3 * (theta - 0.1 * t))
    #     T = np.zeros((Ny, Nx))
    #     ax.set_title(f"t = {t}")
    #     for i in range((Nx) * (Ny)):
    #         l = f.readline().split()
    #         T[int(l[1]), int(l[0])] = float(l[2])

    # X_T, Y_T = np.meshgrid(np.linspace(0, 2/3, endpoint=True, num=(Nx)), np.linspace(0, 1, endpoint=True, num=(Ny)))

    # cs_T = ax.contourf(X_T, Y_T, T, cmap='bwr', levels=levels_Rehw, extend='both')
    # ax.fill(r * np.cos(theta) + 1/3, r * np.sin(theta) + 1/3, color='grey')
    # ax.fill(1/25 * np.cos(theta) + 1/3, 1/25 * np.sin(theta) + 1/3, color='grey')

    # with open("Project/output/U-%d.txt"%(diter*(t+1)), 'r') as f:
    #     l = f.readline().split()
    #     Nx = int(l[2])
    #     Ny = int(l[5])
    #     h = float(l[8])
    #     t = float(l[11])
    #     r = 1/5 * np.cos(3 * (theta - 0.1 * t))
    #     T = np.zeros((Ny, Nx))
    #     ax.set_title(f"t = {t}")
    #     for i in range((Nx-1) * (Ny-1)):
    #         l = f.readline().split()
    #         T[int(l[1]), int(l[0])] = float(l[2])


    # X_T, Y_T = np.meshgrid(np.linspace(0, 2/3, endpoint=True, num=(Nx)), np.linspace(0, 1, endpoint=True, num=(Ny)))

    # cs_T = ax.contourf(X_T, Y_T, T, cmap='cividis', levels=levels_U, extend='min')
    # ax.fill(r * np.cos(theta) + 1/3, r * np.sin(theta) + 1/3, color='grey')
    # ax.fill(1/25 * np.cos(theta) + 1/3, 1/25 * np.sin(theta) + 1/3, color='grey')

    # with open("Project/output/Reh-%d.txt"%(diter*(t+1)), 'r') as f:
    #     l = f.readline().split()
    #     Nx = int(l[2])
    #     Ny = int(l[5])
    #     h = float(l[8])
    #     t = float(l[11])
    #     r = 1/5 * np.cos(3 * (theta - 0.1 * t))
    #     T = np.zeros((Ny, Nx))
    #     ax.set_title(f"t = {t}")
    #     for i in range((Nx) * (Ny)):
    #         l = f.readline().split()
    #         T[int(l[1]), int(l[0])] = float(l[2])

    # X_T, Y_T = np.meshgrid(np.linspace(0, 2/3, endpoint=True, num=(Nx)), np.linspace(0, 1, endpoint=True, num=(Ny)))

    # cs_T = ax.contourf(X_T, Y_T, T, cmap='bwr', levels=levels_Reh, extend='both')
    # ax.fill(r * np.cos(theta) + 1/3, r * np.sin(theta) + 1/3, color='grey')
    # ax.fill(1/25 * np.cos(theta) + 1/3, 1/25 * np.sin(theta) + 1/3, color='grey')
    

anim = FuncAnimation(fig, animate, frames=n_iter, interval=120)

# plt.show()

anim.save('T_350_mix_final.gif')