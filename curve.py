from calculationmethods import *
from scipy.integrate import simps
import numpy as np

sqang = angtocm**2


class Curve:
    def __init__(self, grid, frequency, squared_dipole_moment, name=""):
        """Create a new state object"""
        self.name = name
        self.frequency = np.array(frequency)
        self.grid = np.array(grid)
        self.squared_dipole_moment = squared_dipole_moment
        self.intensity = einstein_coefficient(self.squared_dipole_moment, self.frequency[-1]*cmtos1)
        self.phi = phase_shift_table(self.grid, self.frequency)
        self.u = np.arange(1e4, 6e5, 1e3)
    
    def coefficient_calc(self, T, order=4):
        sigma_s = np.array([simps(self.grid * np.sin(self.phi/v), self.grid) for v in self.u]) * sqang
        sigma_b = np.array([simps(self.grid * (1 - np.cos(self.phi/v)), self.grid) for v in self.u]) * sqang
        weighted_v = np.array([v * maxwell(v, T) for v in self.u])
        k_s = simps(sigma_s * weighted_v, self.u)
        k_b = simps(sigma_b * weighted_v, self.u)
        return k_b + 1j * k_s


    def __str__(self, order=4):
        result = f"States: {self.name}\n"
        result += f"Intensity: {self.intensity}\n"
        result += f"Unperturbed frequency: {self.frequency[-1]}\n"
        return result
