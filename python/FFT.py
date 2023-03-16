import numpy as np
from matplotlib import pyplot as plt
import scipy.fft as fft

def Gaussian(X, sigma):
    return np.exp(-np.power(X/sigma, 2))

def fourierGaussian(j, sigma):
    return np.sqrt(np.pi) * np.exp(-np.power(np.pi * sigma * j, 2))

sigma = 1/16
h = 1/8 * sigma

N = int(1/h)

X = np.linspace(-0.5, 0.5, num=N, endpoint=False)
U = Gaussian(X, sigma)
yf = fft.fft(U)
xf = fft.fftfreq(N, h)
xf = fft.fftshift(xf)
yplot = np.sqrt(np.pi) * sigma * fft.fftshift(yf)
yreal = fourierGaussian(xf, sigma)



plt.figure()
plt.plot(xf, np.abs(yplot))
plt.plot(xf, yreal)
plt.show()
