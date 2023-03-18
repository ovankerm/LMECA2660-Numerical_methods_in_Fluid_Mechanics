import numpy as np
from matplotlib import pyplot as plt

CFLs = np.linspace(1.7, 1.8, num = 100)

kh = np.linspace(0, np.pi, num=100)
k_star_h = -1/6 * (np.exp(-2 * 1j * kh) - 6 * np.exp(-1j * kh) + 2 * np.exp(1j * kh) + 3)

max_CFL = 0

for i, CFL in enumerate(CFLs):
    rho = np.abs(1 + CFL * k_star_h + 1/2 * np.power(CFL * k_star_h, 2) + 1/6 * np.power(CFL * k_star_h, 3) + 1/24 * np.power(CFL*k_star_h, 4))
    if all(ele <= 1 for ele in rho) :
        max_CFL = CFL

print(max_CFL)

# k_star_h *= max_CFL

k_star_h *= 0.5



X, Y = np.meshgrid(np.linspace(-3, 1, endpoint=True, num=1000), np.linspace(-3, 3, endpoint=True, num=1000))

Z = X + 1j * Y

rho = np.abs(1 + Z + 1/2 * np.power(Z, 2) + 1/6 * np.power(Z, 3) + 1/24 * np.power(Z, 4))



fig, ax = plt.subplots(figsize=(5, 7))
ax.grid()
l2, = ax.plot(k_star_h.real, k_star_h.imag, 'k-')
cs = ax.contour(X, Y, rho, colors='k', linestyles='dashed', levels=[1])
h1,l1 = cs.legend_elements()
ax.legend([h1[0], l2], ['RK4', 'line'])
plt.show()