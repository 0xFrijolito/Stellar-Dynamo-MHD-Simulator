import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("output.csv", delimiter=",", skiprows=1)

x = data[:, 0]
y = data[:, 1]
t = data[:, 2]

nx = len(np.unique(x))
ny = len(np.unique(y))
z = t.reshape((ny, nx))

plt.figure(figsize=(6, 4))
plt.imshow(z, origin='lower', cmap='viridis', aspect='auto')
plt.colorbar(label='Temperature')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Heatmap')
plt.show()