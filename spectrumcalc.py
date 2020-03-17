import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad, simps
from calculationmethods import *
from curve import Curve

with open("test_set.txt") as f:
    lines = [line.strip().split() for line in f.readlines()]
    r = np.array([float(line[0]) for line in lines])
    omega = np.array([float(line[1]) for line in lines])
    curve = Curve(r, omega, 1e-36)
    print(curve.omega0)
    rho = 1/np.linspace(5e7, 5e5, 513)
    print(rho)
    x0 = np.tan(np.linspace(-(np.pi-0.00021)/2.01, (np.pi-0.0002)/2.01, 513)) * angtosm
    xx, rr = np.meshgrid(x0, rho)
    s = np.arange(10, 120, 10) * angtosm
    integrand_real = np.zeros((len(x0), len(rho)))
    sums = np.zeros(len(rho))
    for i in range(len(rho)):
        for j in range(len(x0)):
            integrand_real[i, j] = rho[i] * (1-np.cos(curve.eta(rho[i], 1e-7, x0[j], 1e4)))
        sums[i] = simps(integrand_real[i,:], x0)
    print(simps(sums, rho))
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(xx, rr, integrand_real)
    plt.show()
