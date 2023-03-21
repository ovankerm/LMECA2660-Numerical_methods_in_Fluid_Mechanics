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

def exactPacket(time, sigma):
    X = np.linspace(-0.5, 1.5, num=500)
    f = lambda X, t : np.exp(-np.power((X-t)/sigma, 2)) * np.cos(16*np.pi*(X-t))
    U = np.zeros_like(X)
    U += f(X, time)
    U += f(X, time+1)
    U += f(X, time-1)
    return (X, U)
   

# fig, ax = plt.subplots(figsize=(10, 7))

# points = [32, 64, 128]
# h_over_sigma = [1/2, 1/4, 1/8]
# schemes = ["E2", "ED", "E4", "I4"]
# linestyles = ["k-", "k--", "k-.", "k:"]

# errors = np.zeros((4, 3))

# for j,p in enumerate(points):
#     for k,s in enumerate(schemes):
#         with open("output/%s_%dpoints-diagnostics.txt"%(s, p)) as f:
#             lines = f.readlines()
#             for i, l in enumerate(lines):
#                 l = l.split()
#                 if float(l[2]) == 0.5:
#                     errors[k, j] = l[11]
#                     break

# for i in range(4):
#     ax.loglog(h_over_sigma, errors[i], linestyles[i], label="Integration with the %s scheme"%schemes[i])
#     ax.scatter(h_over_sigma, errors[i], c="black", marker='x')

# alphas = np.zeros(4)
# for i in range(4):
#     alphas[i] = (np.log(errors[i,2]/errors[i,0]))/(np.log(h_over_sigma[2]/h_over_sigma[0])) 

# print(alphas)
# ax.legend()
# ax.grid(which="both")
# ax.set_title(r"Evolution of the global error with h/$\sigma$")
# ax.set_xlabel(r"h/$\sigma$")
# ax.set_ylabel("Global error")

# fig.savefig("images/methods_convergence.eps", format='eps')



#----------Plot Solutions----------#

# colors = ["", "black", "blue"]
# points = [32, 64, 128]
# schemes = ["E2", "E4", "ED", "I4"]

# Us = np.empty((3, 4, 2), dtype=np.ndarray)
# Xs = np.empty((3, 4, 2), dtype=np.ndarray)
# exacts = np.empty((2, 2), dtype=np.ndarray)

# for k in range(4):
#     for m in range(3):
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

#                 Us[m, k, j-1] = U
#                 Xs[m, k, j-1] = X_tot

#                 exacts[j-1, 0], exacts[j-1, 1] = exactGaussian(0.5 * j, sigma)
                

# for i in range(3):
#     fig, ax = plt.subplots(2, 2, figsize=(13, 10))
#     fig.suptitle('Simulations using %d points (h/L = %.4e)'%(points[i], 1/points[i]), y=0.93)
#     for j in range(4):
#         for k in range(1, 3):
#             ax[j//2, j%2].set_xlim((-0.5, 1.5))
#             ax[j//2, j%2].set_ylim((-0.7, 1.1))
#             ax[j//2, j%2].grid(visible=True)
#             ax[j//2, j%2].set_title("%s scheme"%(schemes[j]))
#             ax[j//2, j%2].plot(Xs[i, j, k-1], Us[i, j, k-1], label="Simulated solution at ct/L = %.1f"%(0.5*k), linewidth=2, color=colors[k])
#             ax[j//2, j%2].plot(exacts[k-1, 0], exacts[k-1, 1], label="Exact solution at ct/L = %.1f"%(0.5*k), linestyle="dotted", color=colors[k])
#             ax[j//2, j%2].legend(loc='lower left')
#             ax[j//2, j%2].set_ylabel("u/U")
#             ax[j//2, j%2].set_xlabel("x/L")
#     fig.savefig("images/solutions_%d_points.eps"%points[i], format='eps')


#----------Non Uniform Grid----------#

# colors = ["", "black", "blue"]
# points = [32, 64, 128]
# schemes = ["E2", "E4", "ED", "I4"]

# Us = np.empty((4, 2), dtype=np.ndarray)
# Xs = np.empty((4, 2), dtype=np.ndarray)
# Vs = np.empty((4, 2), dtype=np.ndarray)
# ETAs = np.empty((4, 2), dtype=np.ndarray)

# for k in range(4):
#     for j in range(1, 3):
#         with open("output/%s_128points_mapping-%.4f.txt"%(schemes[k], 0.5*j)) as f:
#             l = f.readline().split()
#             sigma = float(l[2])
#             N = int(l[5])
#             h = float(l[8])
                
#             X = np.zeros(N)
#             U = np.zeros(N)
#             ETA = np.zeros(N)
#             V = np.zeros(N)

#             for i in range(N):
#                 l = f.readline().split()
#                 X[i] = float(l[5])
#                 U[i] = float(l[8])
#                 ETA[i] = float(l[11])
#                 V[i] = float(l[14])


#             X_tot = np.append(X, X+1)
#             U = np.append(U, U)
#             ETA_tot = np.append(ETA, ETA+1)
#             V = np.append(V, V)

#             Us[k, j-1] = U
#             Xs[k, j-1] = X_tot
#             Vs[k, j-1] = V
#             ETAs[k, j-1] = ETA_tot
                
# fig, ax = plt.subplots(2, 2, figsize=(13, 10))
# fig.suptitle('Simulations using 128 points (h/L = %.4e) and a non uniform grid'%(1/128), y=0.93)
# for j in range(4):
#     for k in range(1, 3):
#         ax[j//2, j%2].set_xlim((-0.5, 1.5))
#         ax[j//2, j%2].set_ylim((-0.8, 1.7))
#         ax[j//2, j%2].grid(visible=True)
#         ax[j//2, j%2].set_title("%s scheme"%(schemes[j]))
#         ax[j//2, j%2].plot(ETAs[j, k-1], Vs[j, k-1], label="v/U at ct/L = %.1f"%(0.5*k), linewidth=2, color=colors[k])
#         ax[j//2, j%2].plot(Xs[j, k-1], Us[j, k-1], label="u/U at ct/L = %.1f"%(0.5*k), linestyle="dotted", color=colors[k])
#         ax[j//2, j%2].legend(loc='lower left')
#         ax[j//2, j%2].set_ylabel("v/U (u/U)")
#         ax[j//2, j%2].set_xlabel(r"$\xi$/L (x/L)")

# fig.savefig("images/solutions_non_uniform_grid.eps", format='eps')


#----------Wave Packet----------#

colors = ["", "black", "blue"]
points = [32, 64, 128]
schemes = ["E2", "E4", "ED", "I4"]

Us = np.empty((4, 2), dtype=np.ndarray)
Xs = np.empty((4, 2), dtype=np.ndarray)
exacts = np.empty((2, 2), dtype=np.ndarray)

for k in range(4):
    for j in range(1, 3):
        with open("output/%s_128points_wave_packet_mapping-%.4f.txt"%(schemes[k], 0.25*j)) as f:
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

            Us[k, j-1] = U
            Xs[k, j-1] = X_tot
            exacts[j-1, 0], exacts[j-1, 1] = exactPacket(0.25 * j, sigma)
                
fig, ax = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle('Simulations using 128 points (h/L = %.4e) of a wave packet'%(1/128), y=0.93)
for j in range(4):
    for k in range(1, 3):
        ax[j//2, j%2].set_xlim((-0.5, 1.5))
        ax[j//2, j%2].set_ylim((-0.8, 1.8))
        ax[j//2, j%2].grid(visible=True)
        ax[j//2, j%2].set_title("%s scheme"%(schemes[j]))
        ax[j//2, j%2].plot(Xs[j, k-1], Us[j, k-1], label="u/U at ct/L = %.2f"%(0.25*k), linestyle="solid", color=colors[k])
        ax[j//2, j%2].plot(exacts[k-1, 0], exacts[k-1, 1], label="Exact solution at ct/L = %.2f"%(0.25*k), linestyle="dotted", color=colors[k])
        ax[j//2, j%2].legend(loc='upper left')
        ax[j//2, j%2].set_ylabel("u/U")
        ax[j//2, j%2].set_xlabel("x/L")

fig.savefig("images/solutions_wave_packet.eps", format='eps')



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
