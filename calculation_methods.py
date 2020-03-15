import numpy as np
from scipy.optimize import curve_fit

mu = 3.63 # reduced mass for Ar-He
R = 8.31E7 # gas constant
hbar = 1.054e-27 # Dirac constant
c = 2.998E10 # speed of light


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


def potential(r, w0, c1, c2, c3, c4, c5, c6):
    """
    Model of potential functions
    """
    coeffs = [c1, c2, c3, c4, c5, c6] # Yeah, it's very stupid
    return np.array([w0 + sum([coeffs[j] / (r[i] ** (3 * j + 3)) for j in range(6)]) for i in range(len(r))])


def interpolation_coefficients(x, w):
    """
    There used notation of curve_fit method
    Look at scipy.optimize.curve_fit docs
    """
    p, cov = curve_fit(potential, x, w,
                       p0=[w[-1], 3E-9, -4E-31, 1E-53, -3E-76, 1E-99, -2E-123], gtol=1e-18)
    return p[0], p[1:]


def j(n, z, s):
    """
    Calculate integral of function
    (x ^ 2 + s ^ 2) ^ (-n)
    from 0 to z
    """
    if n == 1:
        return z / np.sqrt(s ** 2 + z ** 2)
    elif n == 2:
        return (np.arctan2(s, z) / 2) + s * z / (2 * (s ** 2 + z ** 2))
    else:
        return s * z ** (n - 1) / (n * (s ** 2 + z ** 2) ** (n / 2)) + (n - 1) * j(n - 2, z, s) / n


def phase_from_zero(x, s, v, coeffs):
    """
    Calculate phase integral
    of interpolated frequency function 
    from 0 to x point along trajectory
    on the distance s from fixed point
    """
    return sum([j(3 * i + 1, x, s) * coeffs[i] * (x ** (-3 * i - 2)) for i in range(len(coeffs))]) / v


def phase(s, l, x0, v, coeffs):
    """
    Calculate phase integral
    of interpolated frequency function 
    from x0 point along trajectory
    to the distance l
    on the distance s from fixed point
    """
    return phase_from_zero(x0+l, s, v, coeffs) - phase_from_zero(x0, s, v, coeffs)


def einstein_coefficient(dm2, omega):
    return 4 * dm2 * omega ** 3 / (3 * hbar * c ** 3)