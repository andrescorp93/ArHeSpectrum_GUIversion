from calculationmethods import *
from scipy.integrate import nquad, romberg, quad, simps
from scipy.special import gamma
from numba import jit, njit
from multiprocessing import Pool
from scipy.interpolate import CubicSpline, RectBivariateSpline
import numpy as np

class Curve:
    def __init__(self, grid, frequency, squared_dipole_moment, name=""):
        """Create a new state object"""
        self.name = name
        self.frequency = frequency
        self.grid = grid
        cs = CubicSpline(self.grid, self.frequency, bc_type=((2, 0), (1, 0)))
        self.r = np.arange(2e-8,5e-7,5e-11)
        self.phi = np.zeros(len(self.r))
        for i in range(len(self.r)):
            rho = self.r[i] * np.ones(len(self.r[i:]))
            self.phi[i] = -2 * simps(cs(self.r[i:], 1)*np.sqrt(+self.r[i:]**2-rho**2), self.r[i:])
        self.squared_dipole_moment = squared_dipole_moment
        self.intensity = einstein_coefficient(self.squared_dipole_moment, self.frequency[-1])
    
    def coefficient_calc(self, T, order=4):
        sigma_s = np.array([simps(self.r * np.sin(self.phi/v), self.r) for v in np.arange(1e4, 6e5, 1e3)])
        sigma_b = np.array([simps(self.r * (1 - np.cos(self.phi/v)), self.r) for v in np.arange(1e4, 6e5, 1e3)])
        weighted_v = np.array([v * maxwell(v, T) for v in np.arange(1e4, 6e5, 1e3)])
        k_s = simps(sigma_s * weighted_v, np.arange(1e4, 6e5, 1e3))
        k_b = simps(sigma_b * weighted_v, np.arange(1e4, 6e5, 1e3))
        return k_b + 1j * k_s


    def __str__(self, order=4):
        result = f"States: {self.name}\n"
        result += f"Intensity: {self.intensity}\n"
        result += f"Unperturbed frequency: {self.frequency[-1]}\n"
        return result
