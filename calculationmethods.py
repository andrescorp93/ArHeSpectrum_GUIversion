import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gamma

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


def approx_coeffs(x, y, order=4, n=10):
    """
    There used Method of least squares
    for polynomials
    """
    x_vec = [x ** (-order*i) for i in range(n+1)]
    m = np.array([[np.dot(x_vec[i], x_vec[j]) for j in range(n+1)] for i in range(n+1)])
    b = np.array([np.dot(y, x_vec[i]) for i in range(n+1)])
    return np.linalg.solve(m, b)


def potential(r, c, order=4):
    """
    Model of potential functions
    """
    return np.array([np.dot([x ** (-order*i) for i in range(len(c))], c) for x in r])


def phase_shift(r, v, coeffs):
    """
    Calculate phase integral
    of interpolated frequency function 
    from 0 to x point along trajectory
    on the distance s from fixed point
    """
    k = [np.sqrt(np.pi)*gamma(2*i+(3/2))/gamma(2*i+2) for i in range(len(coeffs))]
    return sum([k[i]*coeffs[i]*(r**(-4*i-3))/v for i in range(len(coeffs))])


def einstein_coefficient(dm2, omega):
    return 4 * dm2 * omega ** 3 / (3 * hbar * c ** 3)
