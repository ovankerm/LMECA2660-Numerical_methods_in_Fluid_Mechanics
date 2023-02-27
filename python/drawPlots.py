import numpy as np
from matplotlib import pyplot as plt

with open("output/test.txt") as f:
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

fig, ax = plt.subplots(figsize=(10, 7))
ax.plot(X, U)
ax.scatter(X, U)
plt.show()
