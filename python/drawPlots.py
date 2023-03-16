import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation 

# def exactGaussian(X, t, sigma):
#     f = lambda X, t : np.exp(-np.power((X-t)/sigma, 2))
#     Xreturn = np.append(X-1, np.append(X, X+1))
#     Ureturn = np.append(f(X, t-1), np.append(Ur, Ur))
#     return Xreturn, Ureturn
   
fig, ax = plt.subplots(figsize=(10, 7))
ax.set_xlim((-0.5, 1.5))
ax.set_ylim((-1, 1.1))
ax.grid()

# lineE2, = ax.plot([], []) 
   
# def initE2(): 
#     return animateE2(0)
   
# def animateE2(i):
#     with open("output/E2-%.4f.txt"%(0.1*i)) as f:
#         l = f.readline().split()
#         L = float(l[2])
#         N = int(l[5])
#         h = float(l[8])
        
#         X = np.zeros(N)
#         U = np.zeros(N)

#         for i in range(N):
#             l = f.readline().split()
#             X[i] = float(l[5])
#             U[i] = float(l[8])

#     lineE2.set_data(X, U)
      
#     return lineE2,
   
# animE2 = FuncAnimation(fig, animateE2, init_func = initE2,
#                      frames = 1000, interval = 100, blit = True)

for j in range(1, 3):
    with open("output/E2-%.4f.txt"%(0.5*j)) as f:
        l = f.readline().split()
        sigma = float(l[2])
        N = int(l[5])
        h = float(l[8])
            
        X = np.zeros(N)
        U = np.zeros(N)

        for i in range(N):
            l = f.readline().split()
            X[i] = float(l[5])
            U[i] = float(l[8])

        X_tot = np.append(X, X+1)
        U = np.append(U, U)

        ax.plot(X_tot, U, label="t = %.4f"%(0.5*j))
        # Xr, Ur = exactGaussian(X, 0.5 * j, sigma)
        # ax.plot(Xr, Ur)

ax.legend()


fig, ax = plt.subplots(figsize=(10, 7))
with open("output/E2-diagnostics.txt") as f:
    lines = f.readlines()
    t = np.zeros(len(lines))
    I = np.zeros(len(lines))
    E = np.zeros(len(lines))
    R = np.zeros(len(lines))
    for i, l in enumerate(lines):
        l = l.split()
        t[i] = l[2]
        I[i] = l[5]
        E[i] = l[8]
        R[i] = l[11]
ax.plot(t, I, label="I")
ax.plot(t, E, label="E")
ax.plot(t, R, label="R")
ax.legend()

plt.show()
