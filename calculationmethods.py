import numpy as np
from scipy.optimize import curve_fit
from scipy.special import gamma

mu = 3.63 # reduced mass for Ar-He
R = 8.31E7 # gas constant
hbar = 1.054e-27 # Dirac constant
c = 2.998E10 # speed of light
cmtos1 = c * 2 * np.pi
angtocm = 1E-8


def power_integral(r, n):
    c = np.sqrt(np.pi)*gamma((n-1)/2)/gamma(n/2)
    return np.power(r, 1-n)*c


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


def approx_coeffs(x, y, k=4, n=10):
    """
    There used Method of least squares
    for polynomials
    """
    if k == 0:
        order = [0, 6, 8, 12]
    else:
        order = [k*i for i in range(m)]
    x_vec = [x ** (-order[i]) for i in range(len(order))]
    m = np.array([[np.dot(x_vec[i], x_vec[j]) for j in range(len(order))] for i in range(len(order))])
    b = np.array([np.dot(y, x_vec[i]) for i in range(len(order))])
    return np.linalg.solve(m, b)


def potential(r, c, k=4, m=10):
    """
    Model of potential functions
    """
    if k == 0:
        order = [0, 6, 8, 12]
    else:
        order = [k*i for i in range(m)]
    x = [r ** (-order[i]) for i in range(len(c))]
    return np.dot(c, x)


def phase_shift(r, coeffs, n=4, m=10):
    """
    Calculate phase integral coefficients
    of interpolated frequency function 
    from 0 to x point along trajectory
    on the distance s from fixed point
    """
    if n == 0:
        order = [0, 6, 8, 12]
    else:
        order = [n*i for i in range(1, m)]
    eta = np.zeros(len(r))
    for k in range(len(coeffs)):
        eta += power_integral(r, order[k]) * coeffs[k]
    return eta * angtocm * cmtos1


def einstein_coefficient(dm2, omega):
    return 4*dm2*omega**3 / (3*hbar*c**3)
