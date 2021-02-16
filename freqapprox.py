import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import curve_fit
from scipy.integrate import simps
from scipy.special import gamma, kn

mu = 3.63  # reduced mass for Ar-He
R = 8.31E7  # gas constant
cmtos1 = 2.998E10 * 2 * np.pi
angtocm = 1E-8


def u(r, u0, A, alpha, C6, C8):
    uex = A * np.exp(-alpha*r)
    u6 = C6 * np.power(r, -6)
    u8 = C8 * np.power(r, -8)
    return u0 + uex + u6 + u8


def approx_coeffs(r, ut):
    popt, pcov = curve_fit(u, r, ut, p0=[ut[-1], -5.6e4, 1.2, 1.9e6, -1.5e7])
    res = {'u0': popt[0], 'A': popt[1], 'alpha': popt[2],
           'C6': popt[3], 'C8': popt[4]}
    return res


def plot_approximation(coeffs):
    r = np.arange(2,10,0.01)
    uapp = u(r, coeffs["u0"], coeffs["A"], coeffs["alpha"], coeffs["C6"], coeffs["C8"])
    plt.plot(r,uapp)

    
def buckingham_phase(r, v, A, alpha):
    A *= cmtos1
    return 2*A*r*kn(1, r)*angtocm / v


def power_integral(r, n):
    k = np.sqrt(np.pi)*gamma((n-1)/2)/gamma(n/2)
    return np.power(r, 1-n)*k


def van_der_waals_phase(r, v, C6, C8):
    C6 *= cmtos1
    C8 *= cmtos1
    phi6 = power_integral(r, 6) * C6
    phi8 = power_integral(r, 8) * C8
    return (phi6+phi8)*angtocm / v


def phase_shift(r, v, up):
    buckinghampart = buckingham_phase(r, v, up["A"], up["alpha"])
    vanderwaalspart = van_der_waals_phase(r, v, up["C6"], up["C8"])
    return buckinghampart + vanderwaalspart


def maxwell(v, T):
    """
    Maxwell distribution
    v - velocity, cm/s
    T - temperature, K
    norm - coefficient of normalization
    """
    norm = ((mu / (2 * np.pi * R * T)) ** 1.5)
    return 4 * np.pi * (v ** 2) * norm * np.exp(-mu * v ** 2 / (2 * R * T))


def calc_coeffs(T, coeffs):
    r = np.arange(1.8, 20, 0.002)
    v = np.arange(1e3, 7e5, 5e2)
    sigmab = np.zeros(len(v))
    sigmas = np.zeros(len(v))
    phaseshift = phase_shift(r, 1000, coeffs) * 1000
    for i in range(len(v)):
        sigmab[i] = simps(r*(1-np.cos(phaseshift/v[i])), r) * angtocm**2
        sigmas[i] = simps(r*(np.sin(phaseshift/v[i])), r) * angtocm**2
    kb = simps(sigmab*v*maxwell(v, T), v)
    ks = simps(sigmas*v*maxwell(v, T), v)
    return kb, ks


with open("results/1s5_2p9.txt", "r") as energy_file:
    energy_txt = [line.strip().split("\t") for line in energy_file.readlines()]
    energies = {energy_txt[0][i]: np.array([float(energy_txt[k][i]) for k in range(2, len(energy_txt))])
                for i in range(len(energy_txt[0]))}
    for k, v in energies.items():
        if k != "R":
            res = approx_coeffs(energies["R"], v)
            print(res)
            plot_approximation(res)
            plt.plot(energies["R"], v)
            plt.show()
            print(calc_coeffs(300, res))
