import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

x = np.loadtxt("sol.txt")
plt.imshow(x, cmap='gray', vmin=0, vmax=1)
#plt.colorbar()
plt.show(block=False)

#To generate a gif save all files
plt.savefig('./images/'+datetime.now().strftime('%H%M%S')+'.png')

plt.pause(1)
plt.close()

