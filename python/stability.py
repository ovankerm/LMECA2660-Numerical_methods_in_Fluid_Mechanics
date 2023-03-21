import numpy as np
from matplotlib import pyplot as plt

CFLs = np.linspace(0, 2 * np.sqrt(2), num = 1000)

kh = np.linspace(0, np.pi, num=100)

k_star_h = np.zeros((4, len(kh)), dtype=np.complex128)

schemes = ["E2", "E4", "ED", "I4"]
k_star_h[0] = np.sin(kh)
k_star_h[1] = 8/6 * np.sin(kh) - 1/6 * np.sin(2 * kh)
k_star_h[2] = -1j/6 * (np.exp(-2 * 1j * kh) - 6 * np.exp(-1j * kh) + 2 * np.exp(1j * kh) + 3)
k_star_h[3] = (3/2 * np.sin(kh))/(1 + 1/2 * np.cos(kh))


max_CFLs = np.zeros(4)

for k in range(4):
    for i, CFL in enumerate(CFLs):
        rho = np.abs(1 - 1j * CFL * k_star_h[k] + 1/2 * np.power(-1j * CFL * k_star_h[k], 2) + 1/6 * np.power(-1j * CFL * k_star_h[k], 3) + 1/24 * np.power(-1j * CFL * k_star_h[k], 4))
        if all(ele <= 1 for ele in rho) :
            max_CFLs[k] = CFL

print(max_CFLs)

k_star_h *= -0.5j

X, Y = np.meshgrid(np.linspace(-3, 1, endpoint=True, num=1000), np.linspace(-3, 3, endpoint=True, num=1000))

Z = X + 1j * Y

rho = np.abs(1 + Z + 1/2 * np.power(Z, 2) + 1/6 * np.power(Z, 3) + 1/24 * np.power(Z, 4))

fig, ax = plt.subplots(figsize=(6, 8))
ax.grid()
lines = np.empty(4, dtype=plt.Line2D)
for i, values in enumerate(k_star_h):
    lines[i], = ax.plot(values.real, values.imag)
cs = ax.contour(X, Y, rho, colors='k', linestyles='dashed', levels=[1])
h1,l1 = cs.legend_elements()
schemes.append('RK4')
lines = np.append(lines, h1[0])
ax.legend(lines, schemes)
plt.show()