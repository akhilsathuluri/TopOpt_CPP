import numpy as np
import matplotlib.pyplot as plt

x = np.loadtxt("sol.txt")
plt.imshow(x, cmap='gray', vmin=0, vmax=1)
plt.show()

