import numpy as np
from matplotlib import pyplot as plt

N = 100
x, dx = np.linspace(0, 2 * np.pi, retstep=True, num=N)
u = np.sin(x)
du = np.cos(x)
dU = np.zeros_like(du)
for i in range(len(dU)):
    dU[i] = 1/(6 * dx) * (u[i-2] - 6 * u[i-1] + 3 * u[i] + 2 * u[(i+1)%N])

# plt.figure()
# plt.plot(x, dU)
# plt.plot(x, du)
# plt.show()

CFLs = np.linspace(0, 3/2 * np.sqrt(2), num = 100)

kh = np.linspace(0, np.pi, num=100)
k_star_h = -1/6 * (np.exp(-2 * 1j * kh) - 6 * np.exp(-1j * kh) + 2 * np.exp(1j * kh) + 3)

max_CFL = 0

for i, CFL in enumerate(CFLs):
    rho = np.abs(1 + CFL * k_star_h + 1/2 * np.power(CFL * k_star_h, 2) + 1/6 * np.power(CFL * k_star_h, 3) + 1/24 * np.power(CFL*k_star_h, 4))
    if all(ele <= 1 for ele in rho) :
        max_CFL = CFL

k_star_h *= max_CFL



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