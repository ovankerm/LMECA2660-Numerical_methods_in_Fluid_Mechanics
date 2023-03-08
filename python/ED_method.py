import numpy as np
from matplotlib import pyplot as plt

N = 100
x,h = np.linspace(0, 2 * np.pi, endpoint=False, retstep=True, num=N)
u = np.sin(x)
du = np.cos(x)
dU = np.zeros_like(u)

for i in range(len(x)):
    dU[i] = 1/(6 * h) * (u[(i-2) % N] - 6 * u[(i-1) % N] + 3 * u[i] + 2 * u[(i+1) % N])


kh = np.linspace(0, np.pi, num=100)
k_star_h = -1j/6 * (np.exp(-2 * 1j * kh) - 6 * np.exp(-1j * kh) + 2 * np.exp(1j * kh) + 3)
real = -1/6 * (np.cos(2 * kh) - 4 * np.cos(kh) + 3)
imag = 1/6 * (np.sin(2 * kh) - 8 * np.sin(kh))

plt.figure()
plt.plot(k_star_h.real, k_star_h.imag)
# plt.plot(x, dU)
plt.show()