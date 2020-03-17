from calculationmethods import *
from scipy.integrate import simps
class Curve:
    def __init__(self, grid, frequency, squared_dipole_moment, name=""):
        """Create a new state object"""
        self.name = name
        self.frequency = frequency
        self.grid = grid
        self.interpolate()
        self.squared_dipole_moment = squared_dipole_moment
        self.intensity = einstein_coefficient(self.squared_dipole_moment, self.omega0)
    
    def interpolate(self):
        self.omega0, self.coefficients = interpolation_coefficients(self.grid, self.frequency)

    def eta(self, rho, v):
        """
        Calculate eta for start position x0
        on the trajectory with length s
        on the distance rho from fixed pont
        with velocity v
        """
        return phase_shift(rho, v, self.coefficients)

    def sigma_calc(self):
        pass

    def coefficient_calc(self):
        pass

    def __str__(self):
        result = f"States: {self.name}\n"
        result += f"Intensity: {self.intensity}\n"
        result += f"Unperturbed frequency: {self.omega0}\n"
        for i in range(len(self.coefficients)):
            result += f"C_{3*(i + 1)} = {self.coefficients[i]}\n"
        return result
