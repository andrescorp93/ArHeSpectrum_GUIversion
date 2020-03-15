import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.integrate import quad, simps
from calculation_methods import *
from curve import Curve
from frequency_calc import angtosm

with open("test_set.txt") as f:
    lines = [line.strip().split() for line in f.readlines()]
    r = np.array([float(line[0]) for line in lines])
    omega = np.array([float(line[1]) for line in lines])
    curve = Curve(r, omega, 1e-36)
    print(curve.omega0)
    rho = np.logspace(0, 5, num=1025, base=2) * angtosm
    print(rho)
