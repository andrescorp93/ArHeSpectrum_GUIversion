from calculationmethods import *
from scipy.integrate import nquad, romberg, quad
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

    def coefficient_calc(self, T):
        phi = lambda r, v: phase_shift(r, v, self.coefficients)
        integrandr = lambda r, v: r*v*maxwell(v, T)*(1-np.cos(phi(r, v)))
        integrandi = lambda r, v: r*v*maxwell(v, T)*np.sin(phi(r, v))
        r = quad(lambda v: romberg(lambda r: integrandr(r, v), 2e-8, 1.5e-7), 1e2, 1e6)[0]
        i = quad(lambda v: romberg(lambda r: integrandi(r, v), 2e-8, 1.5e-7), 1e2, 1e6)[0]
        return r + 1j*i

    def __str__(self):
        result = f"States: {self.name}\n"
        result += f"Intensity: {self.intensity}\n"
        result += f"Unperturbed frequency: {self.omega0}\n"
        for i in range(len(self.coefficients)):
            result += f"C_{3*(i + 2)} = {self.coefficients[i]}\n"
        return result
