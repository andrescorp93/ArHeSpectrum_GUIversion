import numpy as np
import matplotlib.pyplot as plt
from calculationmethods import *
from scipy.interpolate import CubicSpline
from scipy.integrate import simps
from scipy.special import gamma

cmtos1 = 2.998E10 * 2 * np.pi
angtocm = 1E-8


def power_integral(r, n):
    c = np.sqrt(np.pi)*gamma((n-1)/2)/gamma(n/2)
    return np.power(r, 1-n)*c


def approx_coeffs(x, y, k=4, m=10):
    if k == 0:
        order = [0, 6, 8, 12]
    else:
        order = [k*i for i in range(m)]
    x_vec = [x ** (-order[i]) for i in range(len(order))]
    m = np.array([[np.dot(x_vec[i], x_vec[j]) for j in range(len(order))] for i in range(len(order))])
    b = np.array([np.dot(y, x_vec[i]) for i in range(len(order))])
    return np.linalg.solve(m, b)


def potential(r, c, k=4, m=10):
    if k == 0:
        order = [0, 6, 8, 12]
    else:
        order = [k*i for i in range(m)]
    x = [r ** (-order[i]) for i in range(len(c))]
    return np.dot(c, x)


def eta_by_coeffs(r, c, n=4, m=10):
    if n == 0:
        order = [0, 6, 8, 12]
    else:
        order = [n*i for i in range(1, m)]
    eta = np.zeros(len(r))
    for k in range(len(c)):
        eta += power_integral(r, order[k]) * c[k]
    return eta * angtocm * cmtos1


def eta_by_spline(r, spline):
    eta = np.zeros(len(r))
    for k in range(len(r)-1):
        eta[k] = -2 * simps(spline(r[k:], 1)*np.sqrt(r[k:]**2-r[k]**2), r[k:])
    return eta * angtocm * cmtos1


def sin_series(x, n):
    summ = x
    q = x
    for i in range(1, n):
        q *= -np.power(x, 2)/(2*i * (2*i+1))
        summ += q
    return summ



def cos_series(x, n):
    summ = np.power(x, 2) / 2
    q = np.power(x, 2) / 2
    for i in range(1, n):
        q *= -np.power(x, 2)/((2*i+2)*(2*i+1))
        summ += q
    return summ


with open("results/1s5_2p9.txt", "r") as energy_file:
    energy_txt = [line.strip().split("\t") for line in energy_file.readlines()]
    energies = {energy_txt[0][i]: np.array([float(energy_txt[k][i]) for k in range(2, len(energy_txt))])
                for i in range(len(energy_txt[0]))}
    velocities = {}
    with open("velocities.txt", "r") as velocity_file:
        velocity_txt = [line.strip().split("\t") for line in velocity_file.readlines()]
        velocities = {velocity_txt[0][i]: np.array([float(velocity_txt[k][i]) for k in range(2, len(velocity_txt))])
                for i in range(len(velocity_txt[0]))}
    v_r = {velocities["1s5"][k]: velocities["R"][k] for k in range(len(velocities["R"]))}
    for k in energies.keys():
        if k == "R":
            energies["R"] = energies[k][:-1] / angtocm
        else:
            energies[k] = energies[k][:-1] / cmtos1
            coefficients1 = approx_coeffs(energies["R"], energies[k], 6, 4)
            coefficients2 = approx_coeffs(energies["R"], energies[k], 0)
            cs = CubicSpline(energies["R"], energies[k], bc_type=((2, 0), (1, 0)))
            v = np.arange(1e4, 6e5, 1e3)
            x = np.arange(2, 9, 0.05)
            plt.plot(x, eta_by_coeffs(x, coefficients2[1:], 0))
            plt.plot(x, eta_by_coeffs(x, coefficients1[1:], 6, 4))
            plt.plot(x, eta_by_spline(x, cs))
            plt.show()
            #sigma_s = [simps(x * np.sin(eta_by_coeffs(x, coefficients2[1:], 0)/u), x) * angtocm**2 for u in v]
            #sigma_b = [simps(x * (1-np.cos(eta_by_coeffs(x, coefficients2[1:], 0)/u)), x) * angtocm**2 for u in v]
