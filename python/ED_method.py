import numpy as np
from matplotlib import pyplot as plt

kh = np.linspace(0, np.pi, num=100)
k_star_h = 0.5 * 1j/6 * (np.exp(-2 * 1j * kh) - 6 * np.exp(-1j * kh) + 2 * np.exp(1j * kh) + 3)

X, Y = np.meshgrid(np.linspace(-3, 1, endpoint=True, num=1000), np.linspace(-3, 3, endpoint=True, num=1000))

Z = X + 1j * Y

rho = np.abs(1 + Z + 1/2 * np.power(Z, 2) + 1/6 * np.power(Z, 3) + 1/24 * np.power(Z, 4))



fig, ax = plt.subplots()
ax.grid()
l2, = ax.plot(k_star_h.real, k_star_h.imag, 'k-')
cs = ax.contour(X, Y, rho, colors='k', linestyles='dashed', levels=[1])
h1,l1 = cs.legend_elements()
ax.legend([h1[0], l2], ['RK4', 'line'])
plt.show()