import numpy as np
from matplotlib import pyplot as plt


kh_pi = np.linspace(0, 1, num=100)
k_star_h = -1j/6 * (np.exp(-2 * 1j * np.pi * kh_pi) - 6 * np.exp(-1j * np.pi * kh_pi) + 2 * np.exp(1j * np.pi * kh_pi) + 3)



fig, ax = plt.subplots(figsize=(8, 7))
ax.grid()
ax.set_title("Error on the amplitude")
ax.plot(kh_pi, np.abs(k_star_h/np.pi), 'k-', label="ED scheme")
ax.plot(kh_pi, kh_pi, 'k:', label="Ideal behaviour")
ax.set_xlabel(r'$\frac{kh}{\pi}$')
ax.set_ylabel(r'$|\frac{k^*h}{\pi}|$')
ax.legend()
fig.savefig("images/ED_amplitude_error.eps", format='eps')

fig, ax = plt.subplots(figsize=(8, 7))
ax.grid()
ax.set_title("Error on the phase")
ax.plot(kh_pi, np.zeros_like(kh_pi), 'k:', label="Ideal behaviour")
ax.plot(kh_pi, np.angle(k_star_h/np.pi), 'k-', label="ED scheme")
ax.set_xlabel(r'$\frac{kh}{\pi}$')
ax.set_ylabel(r'$arg(\frac{k^*h}{\pi})$')
ax.legend()
fig.savefig("images/ED_phase_error.eps", format='eps')

fig, ax = plt.subplots(figsize=(8, 7))
ax.grid()
ax.set_title("Relative error on the amplitude")
ax.plot(kh_pi, np.abs(k_star_h/(np.pi*kh_pi)), 'k-', label="ED scheme")
ax.plot(kh_pi, np.ones_like(kh_pi), 'k:', label="Ideal behaviour")
ax.set_xlabel(r'$\frac{kh}{\pi}$')
ax.set_ylabel(r'$|\frac{k^*h}{kh}|$')
ax.set_ylim(-0.05, 1.05)
ax.legend()
fig.savefig("images/ED_relative_error.eps", format='eps')

plt.show()