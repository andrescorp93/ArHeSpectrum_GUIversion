import matplotlib.pyplot as plt
import numpy as np
import os
from scipy.integrate import nquad, simps, dblquad
from scipy.optimize import curve_fit

mu = 3.63
R = 8.31E7
n = np.array([2, 3 * np.pi / 8, 32 / 35, 63 * np.pi / 256,
              2048 / 3003, 6435 * np.pi / 32768])
T = np.arange(300, 850, 50)


def potential(r, w0, c1, c2, c3, c4, c5, c6):
    c = [c1, c2, c3, c4, c5, c6]
    return np.array([w0 + sum([c[j] / (r[i] ** (3 * j + 3)) for j in range(6)]) for i in range(len(r))])


def interpolation_coeffs(r, w):
    p, cov = curve_fit(potential, r, w,
                       p0=[w[-1], 3E-9, -4E-31, 1E-53, -3E-76, 1E-99, -2E-123], gtol=1e-18)
    return p[0], p[1:]


def eta(r, v, c):
    return sum([n[i] * c[i] * (r ** (-3 * i - 2)) for i in range(len(c))]) / v


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


os.chdir("results")
with open("1s5_2p10.txt", "r") as freq_file:
    lines = [l.strip().split("\t") for l in freq_file.readlines()]
    frequencies = {"R": np.array([float(lines[i][0]) for i in range(2, len(lines))])}
    for i in range(len(lines[0])):
        frequencies[lines[0][i]] = {"Intensity": float(lines[1][i + 1]),
                                    "Frequency Profile": np.array([float(lines[j][i + 1]) for j in range(2, len(lines))])}

    spectruminfo = {}
    for k in frequencies.keys():
        if k != "R":
            w0, c = interpolation_coeffs(frequencies["R"], frequencies[k]["Frequency Profile"])
            spectruminfo[k] = {"Zero Frequency": w0,
                               "Coefficients": c,
                               "<i|D|f><f|D|i>": frequencies[k]["Intensity"] * 1e-36}
    
    hbar = 1.054e-27
    c = 2.998E10
    aij = 0
    for v in spectruminfo.values():
        aij += 4 * v["<i|D|f><f|D|i>"] * v["Zero Frequency"] ** 3 / (3 * hbar * c ** 3)

    print(aij / 3)
