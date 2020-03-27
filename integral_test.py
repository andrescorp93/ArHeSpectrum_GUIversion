from scipy.integrate import dblquad, nquad, simps, romberg, quad
from scipy.special import gamma
import numpy as np
from calculationmethods import phase_shift, maxwell
import matplotlib.pyplot as plt
import time


p = np.array([1.96e-9, -5.04e-31, 1.41e-53, -1.79e-76, 1.09e-99])
t = 300
rho = np.arange(2e-8, 1e-7, 2e-10)
u = np.arange(1e2, 1e6, 5e4)
phi = lambda r, v: phase_shift(r, v, p)
integrandr = lambda r, v: r*v*maxwell(v, t)*(1-np.cos(phi(r, v)))
integrandi = lambda r, v: r*v*maxwell(v, t)*np.sin(phi(r, v))
for w in u:
    plt.plot(rho, integrandr(rho,w))
plt.show()