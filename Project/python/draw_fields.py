from matplotlib import pyplot as plt
import numpy as np

diter = 1310
n_iter = 1000

data = [10 * diter, 100 * diter, 1000 * diter]
theta = np.linspace(0, 2 * np.pi, num = 100)
levels_omega = np.linspace(-2, 2, num=200)
levels_T = np.linspace(0, 0.01, num=200)
levels_U = np.linspace(0, 0.05, num=200)

mixer = True

size = (15, 8)
ratio = 3/2

fig_T, ax_T = plt.subplots(1, 3, figsize=size)
fig_w, ax_w = plt.subplots(1, 3, figsize=size)
fig_U, ax_U = plt.subplots(1, 3, figsize=size)


for i in range(len(data)):
    with open(f"Project/output/T-{data[i]}.txt", 'r') as f:
        l = f.readline().split()
        Nx = int(l[2])
        Ny = int(l[5])
        h = float(l[8])
        t = float(l[11])
        r = 1/5 * np.cos(3 * (theta - 0.1 * t))
        T = np.zeros((Ny-1, Nx-1))
        for j in range((Nx-1) * (Ny-1)):
            l = f.readline().split()
            T[int(l[1]) - 1, int(l[0]) - 1] = float(l[2])
    X_T, Y_T = np.meshgrid(np.linspace(h/2, 2/3 - h/2, endpoint=True, num=(Nx-1)), np.linspace(h/2, 1 - h/2, endpoint=True, num=(Ny-1)))
    cs_T = ax_T[i].contourf(X_T, Y_T, T, cmap='inferno', levels=levels_T, extend='both')
    if mixer:
        ax_T[i].fill(r * np.cos(theta) + 1/3, r * np.sin(theta) + 1/3, color='grey')
        ax_T[i].fill(1/25 * np.cos(theta) + 1/3, 1/25 * np.sin(theta) + 1/3, color='grey')
    x_left, x_right = ax_T[i].get_xlim()
    y_low, y_high = ax_T[i].get_ylim()
    ax_T[i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax_T[i].set_title(r'$\frac{{U t}}{{H}} = {}$'.format(t))

    with open(f"Project/output/w-{data[i]}.txt", 'r') as f:
        l = f.readline().split()
        Nx = int(l[2])
        Ny = int(l[5])
        h = float(l[8])
        t = float(l[11])
        r = 1/5 * np.cos(3 * (theta - 0.1 * t))
        T = np.zeros((Ny, Nx))
        for j in range((Nx) * (Ny)):
            l = f.readline().split()
            T[int(l[1]), int(l[0])] = float(l[2])
    X_T, Y_T = np.meshgrid(np.linspace(0, 2/3, endpoint=True, num=(Nx)), np.linspace(0, 1, endpoint=True, num=(Ny)))
    cs_w = ax_w[i].contourf(X_T, Y_T, T, cmap='bwr', levels=levels_omega, extend='both')
    if mixer:
        ax_w[i].fill(r * np.cos(theta) + 1/3, r * np.sin(theta) + 1/3, color='grey')
        ax_w[i].fill(1/25 * np.cos(theta) + 1/3, 1/25 * np.sin(theta) + 1/3, color='grey')
    x_left, x_right = ax_w[i].get_xlim()
    y_low, y_high = ax_w[i].get_ylim()
    ax_w[i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax_w[i].set_title(r'$\frac{{U t}}{{H}} = {}$'.format(t))

    with open(f"Project/output/U-{data[i]}.txt", 'r') as f:
        l = f.readline().split()
        Nx = int(l[2])
        Ny = int(l[5])
        h = float(l[8])
        t = float(l[11])
        r = 1/5 * np.cos(3 * (theta - 0.1 * t))
        T = np.zeros((Ny, Nx))
        for j in range((Nx-1) * (Ny-1)):
            l = f.readline().split()
            T[int(l[1]), int(l[0])] = float(l[2])
    X_T, Y_T = np.meshgrid(np.linspace(0, 2/3, endpoint=True, num=(Nx)), np.linspace(0, 1, endpoint=True, num=(Ny)))
    cs_U = ax_U[i].contourf(X_T, Y_T, T, cmap='coolwarm', levels=levels_U, extend='min')
    if mixer:
        ax_U[i].fill(r * np.cos(theta) + 1/3, r * np.sin(theta) + 1/3, color='grey')
        ax_U[i].fill(1/25 * np.cos(theta) + 1/3, 1/25 * np.sin(theta) + 1/3, color='grey')
    x_left, x_right = ax_U[i].get_xlim()
    y_low, y_high = ax_U[i].get_ylim()
    ax_U[i].set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)
    ax_U[i].set_title(r'$\frac{{U t}}{{H}} = {}$'.format(t))


fig_T.subplots_adjust(right=0.8)
cbar_ax = fig_T.add_axes([0.85, 0.25, 0.05, 0.5])
b = fig_T.colorbar(cs_T, cax=cbar_ax)

fig_w.subplots_adjust(right=0.8)
cbar_ax = fig_w.add_axes([0.85, 0.25, 0.05, 0.5])
b = fig_w.colorbar(cs_w, cax=cbar_ax)

fig_U.subplots_adjust(right=0.8)
cbar_ax = fig_U.add_axes([0.85, 0.25, 0.05, 0.5])
b = fig_U.colorbar(cs_U, cax=cbar_ax)

fig_T.savefig("Project/images/T_mix.png", format='png', dpi=300)
fig_w.savefig("Project/images/w_mix.png", format='png', dpi=300)
fig_U.savefig("Project/images/U_mix.png", format='png', dpi=300)

plt.show()