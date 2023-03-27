import numpy as np
from matplotlib import pyplot as plt
import scipy.fft as fft

epsilon_machine = np.finfo(np.float64).eps

def Gaussian(X, sigma):
    return np.exp(-np.power(X/sigma, 2))

def fourierGaussian(j, sigma):
    U = np.sqrt(np.pi) * sigma * np.exp(-np.power(np.pi * sigma * j, 2))
    U = np.where(U < epsilon_machine, epsilon_machine, U)
    return U

fig, ax = plt.subplots(1, 2, figsize=(13, 7))

for i in range(2):
    sigma = 1/(4**(i+1))
    h = 1/8 * sigma

    N = int(1/h)

    X = np.linspace(-0.5, 0.5, num=N, endpoint=False)
    U = Gaussian(X, sigma)
    U_hat = 1/N * fft.fft(U)
    j = fft.fftfreq(N, h)
    j = fft.fftshift(j)
    U_hat = fft.fftshift(U_hat)
    U_real = fourierGaussian(j, sigma)

    ax[i].plot(j, np.log(np.abs(U_hat)), 'xk-', label='Discrete Fourier transform in a periodic domain')
    ax[i].plot(j, np.log(U_real), 'k--', label='Fourier transform in an unbounded domain')
    ax[i].hlines(np.log(epsilon_machine), j[0], j[-1], colors='red', label='Machine epsilon')
    ax[i].grid()
    ax[i].legend(loc='lower left')
    ax[i].set_xlabel("j")
    ax[i].set_title(r'Fourier transform using $\frac{L}{\sigma}$ = %d'%(4**(i+1)))
    ax[i].set_ylabel(r'$log(\hat{F}_j)$')
ax[0].set_ylim(ax[1].get_ylim())
fig.savefig("images/FFT.eps", format='eps')
plt.show()
