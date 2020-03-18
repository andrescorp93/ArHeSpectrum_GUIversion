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

    def interpolated(self):
        c1, c2, c3, c4, c5, c6 = self.coefficients
        return potential(self.grid, self.omega0, c1, c2, c3, c4, c5, c6)

    def eta(self, rho, v):
        """
        Calculate eta
        on the distance rho from fixed pont
        with velocity v
        """
        return phase_shift(rho, v, self.coefficients)

    def sigma_calc(self, v):
        rho = np.arange(2e-8, 1e-7, 1e-10)
        phi = self.eta(rho, v)
        return simps(rho * (np.ones(len(rho))-np.cos(phi)+1j*np.sin(phi)), rho)

    def coefficient_calc(self):
        pass

    def __str__(self):
        result = f"States: {self.name}\n"
        result += f"Intensity: {self.intensity}\n"
        result += f"Unperturbed frequency: {self.omega0}\n"
        for i in range(len(self.coefficients)):
            result += f"C_{3*(i + 1)} = {self.coefficients[i]}\n"
        return result
