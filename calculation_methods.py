import numpy as np
from scipy.optimize import curve_fit

mu = 3.63
R = 8.31E7
n = np.array([2, 3 * np.pi / 8, 32 / 35, 63 * np.pi / 256,
              2048 / 3003, 6435 * np.pi / 32768])
hbar = 1.054e-27
c = 2.998E10


def maxwell(v, T):
    norm = ((mu / (2 * np.pi * R * T)) ** 1.5)
    return 4 * np.pi * (v ** 2) * norm * np.exp(-mu * v ** 2 / (2 * R * T))


def lorentz(w, A, s, b, h = 0):
    """
    Функция Лоренца для генерации спектров
    """

    return b * A / (2 * np.pi * ((w - s) ** 2 + (b / 2) ** 2)) + h


def linear_func(x, k, b=0):
    return np.array([k * x[i] + b for i in range(len(x))])


def potential(r, w0, c1, c2, c3, c4, c5, c6):
    c = [c1, c2, c3, c4, c5, c6]
    return np.array([w0 + sum([c[j] / (r[i] ** (3 * j + 3)) for j in range(6)]) for i in range(len(r))])


def interpolation_coefficients(x, w):
    p, cov = curve_fit(potential, x, w,
                       p0=[w[-1], 3E-9, -4E-31, 1E-53, -3E-76, 1E-99, -2E-123], gtol=1e-18)
    return p[0], p[1:]


def j(n, z, s):
    if n == 1:
        return z / np.sqrt(s ** 2 + z ** 2)
    elif n == 2:
        return (np.arctan2(s, z) / 2) + s * z / (2 * (s ** 2 + z ** 2))
    else:
        return s * z ** (n - 1) / (n * (s ** 2 + z ** 2) ** (n / 2)) + (n - 1) * j(n - 2, z, s) / n


def eta(x, s, v, c):
    return 2 * sum([j(3 * i + 1, x, s) * c[i] * (x ** (-3 * i - 2)) for i in range(len(c))]) / v


def einstein_coefficient(dm2, omega):
    return 4 * dm2 * omega ** 3 / (3 * hbar * c ** 3)