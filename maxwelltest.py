from calculationmethods import maxwell
import numpy as np
import matplotlib.pyplot as plt

v = np.arange(1e4, 6e5, 1e3)
plt.plot(v, maxwell(v, 300))
plt.plot(v, maxwell(v, 500))
plt.plot(v, maxwell(v, 700))
plt.plot(v, maxwell(v, 900))
plt.plot(v, maxwell(v, 1100))
plt.show()
