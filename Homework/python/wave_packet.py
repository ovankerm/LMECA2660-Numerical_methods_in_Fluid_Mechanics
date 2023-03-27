import numpy as np
import scipy.fft as fft
from matplotlib import pyplot as plt

def wavePacket(X, sigma):
    return np.exp(-np.power(X/sigma, 2)) * np.cos(16 * np.pi * X)

sigma = 1/16
h = 1/8 * sigma
N = int(1/h)

X = np.linspace(-0.5, 0.5, num=N, endpoint=False)
U = wavePacket(X, sigma)
U_hat = 1/N * fft.fft(U)
j = fft.fftfreq(N, h)
j = fft.fftshift(j)
U_hat = fft.fftshift(U_hat)

schemes = ["E2", "E4", "ED", "I4"]
linestyles = ["k-", "k--", "k-.", "k:"]

kh = np.linspace(0, np.pi, num=(len(j)//2))

cgs = np.empty((4, len(kh)), dtype=np.complex128)

cgs[0, :] = np.cos(kh)
cgs[1, :] = 8/6 * np.cos(kh) - 1/3 * np.cos(2 * kh)
cgs[2, :] = 1/3 * (-np.exp(-2j * kh) + 3 * np.exp(-1j * kh) + np.exp(1j * kh))
cgs[3, :] = 3/2 * (np.cos(kh) + 1/2)/(np.power(1 + 1/2 * np.cos(kh), 2))

fig, ax = plt.subplots(figsize=(7, 5))
ax.set_title("FFT coefficients for wave packet")
ax.plot(j, np.abs(U_hat), "k-", label="FFT")
ax.legend()
ax.grid()
ax.set_xlabel('j')
ax.set_ylabel(r'|$\hat{f}(x)$|')
fig.savefig("images/FFT_wave_packet.eps", format='eps')

fig, ax = plt.subplots(figsize=(7, 5))
ax.set_title("Group velocities of each schemes")
for i in range(4):
    ax.plot(j[len(j)//2:], cgs[i].real, linestyles[i], label=schemes[i])
ax.legend()
ax.grid()
ax.set_xlabel('j')
ax.set_ylabel(r'$\frac{c^*_g}{c}$', fontsize=15)
fig.savefig("images/group_velocities.eps", format='eps')

plt.show()