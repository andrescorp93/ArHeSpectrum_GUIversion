import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.integrate import quad, simps


def potential(x, w0, c1, c2, c3, c4, c5, c6):
    c = [c1, c2, c3, c4, c5, c6]
    return np.array([w0 + sum([c[j] / (x[i] ** (3 * j + 3)) for j in range(6)]) for i in range(len(x))])


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
