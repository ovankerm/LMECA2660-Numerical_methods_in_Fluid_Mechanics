import numpy as np
from matplotlib import pyplot as plt


kh_pi = np.linspace(0, 1, num=100)
k_star_h = -1j/6 * (np.exp(-2 * 1j * np.pi * kh_pi) - 6 * np.exp(-1j * np.pi * kh_pi) + 2 * np.exp(1j * np.pi * kh_pi) + 3)



fig, ax = plt.subplots()
ax.grid()
ax.plot(kh_pi, np.abs(k_star_h/np.pi), 'k-', label="amplitude")
ax.plot(kh_pi, kh_pi, 'k-.', label="ideal")
# ax.plot(kh_pi, np.angle(k_star_h/np.pi), 'k--', label="phase")
ax.legend()
plt.show()