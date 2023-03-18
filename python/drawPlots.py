import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation 

def exactGaussian(time, sigma):
    X = np.linspace(-0.5, 1.5, num=500)
    f = lambda X, t : np.exp(-np.power((X-t)/sigma, 2))
    U = np.zeros_like(X)
    U += f(X, time)
    U += f(X, time+1)
    U += f(X, time-1)
    return (X, U)
   

fig, ax = plt.subplots(figsize=(10, 7))

points = [32, 64, 128]
h_over_sigma = [1/2, 1/4, 1/8]
schemes = ["E2", "ED", "E4", "I4"]
linestyles = ["k-", "k--", "k-.", "k:"]

errors = np.zeros((4, 3))

for j,p in enumerate(points):
    for k,s in enumerate(schemes):
        with open("output/%s_%dpoints-diagnostics.txt"%(s, p)) as f:
            lines = f.readlines()
            for i, l in enumerate(lines):
                l = l.split()
                if float(l[2]) == 0.5:
                    errors[k, j] = l[11]
                    break

for i in range(4):
    ax.loglog(h_over_sigma, errors[i], linestyles[i], label="Integration with the %s scheme"%schemes[i])
    ax.scatter(h_over_sigma, errors[i], c="black", marker='x')

ax.legend()
ax.grid(which="both")
ax.set_title("Order of convergence of each integration scheme")
ax.set_xlabel(r"$\frac{h}{\sigma}$")
ax.set_ylabel("Global error")



#----------Plot Solutions----------#

# colors = ["", "black", "blue"]
# points = [32, 64, 128]
# schemes = ["E2", "E4", "ED", "I4"]

# for k in range(4):
#     for m in range(3):
#         fig, ax = plt.subplots(figsize=(10, 7))
#         ax.set_xlim((-0.5, 1.5))
#         ax.set_ylim((-1, 1.1))
#         ax.grid()
#         ax.set_title("Solution for the %s scheme and %d points"%(schemes[k], points[m]))
#         for j in range(1, 3):
#             with open("output/%s_%dpoints-%.4f.txt"%(schemes[k], points[m], 0.5*j)) as f:
#                 l = f.readline().split()
#                 sigma = float(l[2])
#                 N = int(l[5])
#                 h = float(l[8])
                    
#                 X = np.zeros(N)
#                 U = np.zeros(N)

#                 for i in range(N):
#                     l = f.readline().split()
#                     X[i] = float(l[5])
#                     U[i] = float(l[8])

#                 X_tot = np.append(X, X+1)
#                 U = np.append(U, U)

#                 ax.plot(X_tot, U, label="simulation ct/L = %.1f"%(0.5*j), linewidth=2, color=colors[j])
#                 Xr, Ur = exactGaussian(0.5 * j, sigma)
#                 ax.plot(Xr, Ur, label="real solution ct/L = %.1f"%(0.5*j), linestyle="dotted", color=colors[j])
                
#         ax.set_xlabel("x/L")
#         ax.set_ylabel("u/U")
#         ax.legend()


#----------Plot Diagnostics----------#

# fig, ax = plt.subplots(3, 1, sharex=True, figsize=(10, 7))
# with open("output/%s_%dpoints-diagnostics.txt"%(schemes[0], points[0])) as f:
#     lines = f.readlines()
#     t = np.zeros(len(lines))
#     I = np.zeros(len(lines))
#     E = np.zeros(len(lines))
#     R = np.zeros(len(lines))
#     for i, l in enumerate(lines):
#         l = l.split()
#         t[i] = l[2]
#         I[i] = l[5]
#         E[i] = l[8]
#         R[i] = l[11]

# labels = ["I", "E", "R"]
# data = [I, E, R]

# for i in range(3):
#     ax[i].plot(t, data[i])
#     ax[i].set_title(labels[i])


# ax[2].set_xlabel("ct/L")
plt.show()







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
