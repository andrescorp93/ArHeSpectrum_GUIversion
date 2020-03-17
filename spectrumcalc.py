import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad, simps
from calculationmethods import *
from curve import Curve
import time

with open("test_set.txt") as f:
    lines = [line.strip().split() for line in f.readlines()]
    r = np.array([float(line[0]) for line in lines])
    omega = np.array([float(line[1]) for line in lines])
    curve = Curve(r, omega, 1e-36)
    rho = np.linspace(2e-8, 1.2e-7, 200)
    print(rho)
    x0 = np.linspace(-1e-6, 1e-6, 400) * angtosm
    xx, rr = np.meshgrid(x0, rho)
    ss = np.arange(1e-7, 1.6e-6, 1e-7)
    gs = []
    for s in ss:
        x = time.time()
        gs.append(curve.g(rho, x0, s, 1e4))
        print(f"{time.time() - x} s")
    plt.plot(ss, np.real(gs))
    plt.plot(ss, np.imag(gs))
    plt.show()

