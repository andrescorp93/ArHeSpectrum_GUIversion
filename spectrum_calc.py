import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import quad, simps
from calculation_methods import *
from curve import Curve

with open("test_set.txt") as f:
    lines = [line.strip().split() for line in f.readlines()]
    r = np.array([float(line[0]) for line in lines])
    omega = np.array([float(line[1]) for line in lines])
    curve = Curve(r, omega, 1e-36)
    print(curve.omega0)
    rho = np.logspace(0, 5, num=1025, base=2) * angtosm
    x0 = np.tan(np.linspace(-(np.pi+0.0001)/4, np.pi/4, 513)) * angtosm
    s = np.arange(8000, 16000, 1000) * angtosm
    for d in s:
        print(d, simps(np.array([simps(rho - rho * np.exp(-1j * curve.eta(rho, d, x, 1e4)), rho) for x in x0]), x0))
        
