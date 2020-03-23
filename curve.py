from calculationmethods import *
from scipy.integrate import simps
class Curve:
    def __init__(self, grid, frequency, squared_dipole_moment, name=""):
        """Create a new state object"""
        self.name = name
        self.omega0 = frequency[-1]
        self.frequency = frequency - self.omega0 * np.ones(len(frequency))
        self.grid = grid
        self.interpolate()
        self.squared_dipole_moment = squared_dipole_moment
        self.intensity = einstein_coefficient(self.squared_dipole_moment, self.omega0)
    
    def interpolate(self):
        self.coefficients = interpolation_coefficients(self.grid, self.frequency)

    def interpolated(self):
        c1, c2, c3, c4, c5, c6 = self.coefficients
        return potential(self.grid, c1, c2, c3, c4, c5, c6)

    def eta(self, rho, v):
        """
        Calculate eta
        on the distance rho from fixed pont
        with velocity v
        """
        return phase_shift(rho, v, self.coefficients)

    def sub_integrand(self, r, T):
        v = np.arange(10, 1.05e6, 5e2)
        phi = self.eta(r, v)
        return simps([maxwell(u, T)*u*r*(1-np.cos(phi)+1j*np.sin(phi)) for u in v], v)

    def coefficient_calc(self, T):
        rho = np.arange(1.9e-8, 1e-7, 2e-9)
        integrand = np.array([self.sub_integrand(r, T) for r in rho])
        return simps(integrand, rho)

    def __str__(self):
        result = f"States: {self.name}\n"
        result += f"Intensity: {self.intensity}\n"
        result += f"Unperturbed frequency: {self.omega0}\n"
        for i in range(len(self.coefficients)):
            result += f"C_{3*(i + 1)} = {self.coefficients[i]}\n"
        return result
