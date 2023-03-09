import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation 
   
fig, ax = plt.subplots(figsize=(10, 7))
ax.set_xlim((-8, 8))
ax.set_ylim((-1, 1.1))
ax.grid()

lineE2, = ax.plot([], []) 
   
def initE2(): 
    return animateE2(0)
   
def animateE2(i):
    with open("output/E2-%.4f.txt"%(0.1*i)) as f:
        l = f.readline().split()
        L = float(l[2])
        N = int(l[5])
        h = float(l[8])
        
        X = np.zeros(N)
        U = np.zeros(N)

        for i in range(N):
            l = f.readline().split()
            X[i] = float(l[5])
            U[i] = float(l[8])

    lineE2.set_data(X, U)
      
    return lineE2,
   
animE2 = FuncAnimation(fig, animateE2, init_func = initE2,
                     frames = 1000, interval = 100, blit = True)


plt.show()
