import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.integrate import quad, simps
from calculation_methods import *

with open("test_set.txt") as f:
    lines = [line.strip().split() for line in f.readlines()]
    r = np.array([float(line[0]) for line in lines])
    omega = np.array([float(line[1]) for line in lines])
    w0, b = interpolation_coefficients(r, omega)
    print(w0)
    a = np.arange(8e-7, 8e-5, 8e-7)
    sigma_real = np.zeros(len(a))
    sigma_imag = np.zeros(len(a))
    z = np.arange(1e-10, 1e-5, 1e-9)
    u = 10000
    for i in range(len(a)):
        sigma_imag[i] = simps(np.array([rho * np.sin(eta(rho, a[i], u, b)) for rho in z]), z)
        sigma_real[i] = simps(np.array([rho * (1 - np.cos(eta(rho, a[i], u, b))) for rho in z]), z)
    plt.plot(sigma_real, sigma_imag)
    plt.show()
