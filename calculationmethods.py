import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gamma, hyp2f1

mu = 3.63 # reduced mass for Ar-He
R = 8.31E7 # gas constant
hbar = 1.054e-27 # Dirac constant
c = 2.998E10 # speed of light
angtosm = 1e-8


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


def potential(r, c1, c2, c3, c4, c5, c6):
    """
    Model of potential functions
    """
    coeffs = [c1, c2, c3, c4, c5, c6] # Yeah, it's very stupid
    return np.array([sum([coeffs[j] / (r[i] ** (3 * j + 3)) for j in range(6)]) for i in range(len(r))])


def interpolation_coefficients(x, w):
    """
    There used notation of curve_fit method
    Look at scipy.optimize.curve_fit docs
    """
    p, cov = curve_fit(potential, x, w,
                       p0=[7.7E-9, -8.5E-31, 2.3E-53, -3.1E-76, 2.0E-99, -4.9E-123], gtol=1e-18)
    return p


def phase_shift(r, v, coeffs):
    """
    Calculate phase integral
    of interpolated frequency function 
    from 0 to x point along trajectory
    on the distance s from fixed point
    """
    k = [np.sqrt(np.pi)*gamma((3*i+2)/2)/gamma((3*i+3)/2) for i in range(len(coeffs))]
    return sum([k[i]*coeffs[i]*(r**(-3*i-2))/v for i in range(len(coeffs))])


def einstein_coefficient(dm2, omega):
    return 4 * dm2 * omega ** 3 / (3 * hbar * c ** 3)
