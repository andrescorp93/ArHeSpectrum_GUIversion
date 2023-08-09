import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gamma
from scipy.integrate import simps

mu = 3.63 # reduced mass for Ar-He
R = 8.31E7 # gas constant
hbar = 1.054e-27 # Dirac constant
c = 2.998E10 # speed of light
cmtos1 = c * 2 * np.pi
angtocm = 1E-8


def maxwell(v, T):
    """
    Maxwell distribution
    v - velocity, cm/s
    T - temperature, K
    norm - coefficient of normalization
    """
    norm = ((mu / (2 * np.pi * R * T)) ** 1.5)
    return 4 * np.pi * (v ** 2) * norm * np.exp(-mu * v ** 2 / (2 * R * T))


def lorentz(w, A, s, b, h = 0):
    """
    Lorentz spectrum function
    """
    return b * A / (2 * np.pi * ((w - s) ** 2 + (b / 2) ** 2)) + h


def linear_func(x, k, b=0):
    return np.array([k * x[i] + b for i in range(len(x))])

def phase_shift_table(r, f):
    """
    Calculate phase integral coefficients
    of tabulated frequency function 
    from rmin to rmax point along trajectory
    on the distance s from fixed point
    """
    eta = np.array([simps(np.array([r[j] * (f[j]-f[-1]) / np.sqrt(r[j]**2 - r[i]**2) if r[j] > r[i] else 0 for j in range(len(r))]), r) for i in range(len(r))])
    return eta * angtocm * cmtos1


def einstein_coefficient(dm2, omega):
    return 4*dm2*omega**3 / (3*hbar*c**3)
