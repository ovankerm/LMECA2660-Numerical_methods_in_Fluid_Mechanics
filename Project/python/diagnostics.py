import numpy as np
from matplotlib import pyplot as plt

N = 1000

dt = 3.809160305e-3

Tavg = np.zeros(N)
sigT = np.zeros(N)
iter = np.zeros(N)
qe = np.zeros(N)
Ts = np.zeros(N)
Rehw = np.zeros(N)
Reh = np.zeros(N)

with open("Project/output/diagnostics.txt", "r") as f:
    for i in range(N):
        l = f.readline().split()
        Tavg[i] = float(l[5])
        sigT[i] = float(l[8])
        iter[i] = float(l[2])
        qe[i] = float(l[11])
        Ts[i] = float(l[14])
        Rehw[i] = float(l[17])
        Reh[i] = float(l[20])

iter *= dt

f = np.poly1d(np.polyfit(iter, qe, 5))

plots = [Tavg, sigT, qe, Ts]
labels = [r"$\frac{(\langle T\rangle (t') - T_{\infty})}{\Delta T}$", r"$\frac{( T_{rms}(t') - T_{\infty})}{\Delta T}$", r"$\frac{\langle q_e\rangle(t')}{q_w}$", r"$\frac{(\langle T\rangle|_{cyl} (t') - T_{\infty})}{\Delta T}$"]

fig, ax = plt.subplots(4, 1, figsize = (12, 9), sharex='all')
for i in range(4):
    ax[i].plot(iter, plots[i], 'k-')
    ax[i].grid()
    ax[i].set_ylabel(labels[i], fontsize=12)
fig.suptitle('Diagnostics', y=0.92)
ax[-1].set_xlabel(r"$t' = \frac{tU}{H}$", fontsize=12)
ax[2].plot(iter, f(iter), 'r', label='Polynomial fit', linewidth=1)
ax[2].legend()

fig.savefig("Project/images/diag_mix.eps", format='eps')

# plt.figure()
# plt.plot(iter, Ts, label='Ts')
# plt.legend()

plots = [Rehw, Reh]
labels = [r'$Re_{h,\omega}$', r'$Re_h$']
fig, ax = plt.subplots(2, 1, figsize = (12, 7), sharex='all')
for i in range(2):
    ax[i].plot(iter, plots[i], 'k-', label=labels[i])
    ax[i].grid()
    ax[i].set_ylabel(labels[i], fontsize=12)
fig.suptitle('Reynolds numbers', y=0.92)
ax[-1].set_xlabel(r"$t' = \frac{tU}{H}$", fontsize=12)

# fig.savefig("Project/images/Reynolds_mix.eps", format='eps')




# plt.show()