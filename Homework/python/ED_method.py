import numpy as np
from matplotlib import pyplot as plt


kh_pi = np.linspace(0, 1, num=100)

k_star_h = -1j/6 * (np.exp(-2 * 1j * np.pi * kh_pi) - 6 * np.exp(-1j * np.pi * kh_pi) + 2 * np.exp(1j * np.pi * kh_pi) + 3)


fig, ax = plt.subplots(figsize=(8, 7))
ax.grid()
ax.set_title(r"Real part of $\frac{k^*h}{\pi}$")
ax.plot(kh_pi, (k_star_h/np.pi).real,'k-', label="For the ED scheme")
ax.plot(kh_pi, kh_pi, 'k:', label="Ideal behaviour")
ax.set_xlabel(r'$\frac{kh}{\pi}$', fontsize=15)
ax.set_ylabel(r'$Re\{\frac{k^*h}{\pi}\}$', fontsize=12)
ax.legend()
fig.savefig("images/ED_real_part.eps", format='eps')

fig, ax = plt.subplots(figsize=(8, 7))
ax.grid()
ax.set_title(r"Imaginary part of $\frac{k^*h}{\pi}$")
ax.plot(kh_pi, (k_star_h/np.pi).imag,'k-', label="For the ED scheme")
ax.plot(kh_pi, np.zeros_like(kh_pi), 'k:', label="Ideal behaviour")
ax.set_xlabel(r'$\frac{kh}{\pi}$', fontsize=15)
ax.set_ylabel(r'$Im\{\frac{k^*h}{\pi}\}$', fontsize=12)
ax.legend()
fig.savefig("images/ED_imag_part.eps", format='eps')

# fig, ax = plt.subplots(figsize=(8, 7))
# ax.grid()
# ax.set_title("Relative error on the amplitude")
# ax.plot(kh_pi, np.abs(k_star_h/(np.pi*kh_pi)), 'k-', label="ED scheme")
# ax.plot(kh_pi, np.ones_like(kh_pi), 'k:', label="Ideal behaviour")
# ax.set_xlabel(r'$\frac{kh}{\pi}$', fontsize=15)
# ax.set_ylabel(r'$|\frac{k^*h}{kh}|$', fontsize=15)
# ax.set_ylim(-0.05, 1.05)
# ax.legend()
# fig.savefig("images/ED_relative_error.eps", format='eps')

plt.show()