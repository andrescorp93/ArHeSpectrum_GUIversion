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
    rho = np.logspace(1, 4, num=513, base=2) * angtosm
    x0 = np.tan(np.linspace(-(np.pi-0.00021)/2, (np.pi-0.0002)/2, 513)) * angtosm
    s = np.arange(0.1, 10, 0.1) * angtosm
    g_re = np.zeros(len(s))
    g_im = np.zeros(len(s))
    for i in range(len(s)):
        g_re[i] = simps(np.array([simps(rho - rho*np.cos(curve.eta(rho, s[i], x, 1e4)), rho) for x in x0]), x0)
        g_im[i] = simps(np.array([simps(rho * np.sin(curve.eta(rho, s[i], x, 1e4)), rho) for x in x0]), x0)
        print(g_re[i], g_im[i])
    plt.plot(s, g_re)
    plt.plot(s, g_im)
    plt.show()
        
