import matplotlib.pyplot as plt
import os
from scipy.integrate import nquad, simps, dblquad
from calculation_methods import *
import curve

T = np.arange(300, 850, 50)

os.chdir("results")
with open("1s5_2p10.txt", "r") as freq_file:
    lines = [l.strip().split("\t") for l in freq_file.readlines()]
    r = np.array([float(lines[i][0]) for i in range(2, len(lines))])
    curves = []
    for i in range(len(lines[0])):
        omega = np.array([float(lines[j][i + 1]) for j in range(2, len(lines))])
        curves.append(curve.Curve(r, omega, float(lines[1][i + 1]) * 1e-36))

    aij = 0
    for curve in curves:
        aij += curve.intensity

    print(aij / 3)
