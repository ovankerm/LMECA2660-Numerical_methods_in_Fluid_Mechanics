import numpy as np
from matplotlib import pyplot as plt

with open("output/E2.txt") as f:
    line = f.readline().split()
    L = float(line[2])
    N = int(line[5])
    h = float(line[8])
    
    X = np.zeros(N)
    U = np.zeros(N)

    for i in range(N):
        line = f.readline().split()
        X[i] = float(line[5])
        U[i] = float(line[8])

exact = np.exp(-(np.power(X - 1, 2)))
fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(X, U)
ax.plot(X, exact)
ax.scatter(X, U)

with open("output/E4.txt") as f:
    line = f.readline().split()
    L = float(line[2])
    N = int(line[5])
    h = float(line[8])
    
    X = np.zeros(N)
    U = np.zeros(N)

    for i in range(N):
        line = f.readline().split()
        X[i] = float(line[5])
        U[i] = float(line[8])

exact = np.exp(-(np.power(X - 1, 2)))
fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(X, U)
ax.plot(X, exact)
ax.scatter(X, U)


with open("output/I4.txt") as f:
    line = f.readline().split()
    L = float(line[2])
    N = int(line[5])
    h = float(line[8])
    
    X = np.zeros(N)
    U = np.zeros(N)

    for i in range(N):
        line = f.readline().split()
        X[i] = float(line[5])
        U[i] = float(line[8])

exact = np.exp(-(np.power(X - 1, 2)))
fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(X, U)
ax.plot(X, exact)
ax.scatter(X, U)
plt.show()
