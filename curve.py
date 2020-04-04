from calculationmethods import *
from scipy.integrate import nquad, romberg, quad
class Curve:
    def __init__(self, grid, frequency, squared_dipole_moment, name=""):
        """Create a new state object"""
        self.name = name
        self.frequency = frequency
        self.grid = grid
        self.interpolate()
        self.squared_dipole_moment = squared_dipole_moment
        self.intensity = einstein_coefficient(self.squared_dipole_moment, self.coefficients[0])
    
    def interpolate(self):
        preres = approx_coeffs(self.grid/angtosm, self.frequency/(2*np.pi*c))
        conv = [preres[i] * (2*np.pi*c) * angtosm**(6*i) for i in range(len(preres))]
        self.coefficients = np.array(conv)

    def interpolated(self):
        return potential(self.grid, self.coefficients)

    def coefficient_calc(self, T):
        phi = lambda r, v: phase_shift(r, v, self.coefficients[1:])
        integrandr = lambda r, v: r*v*maxwell(v, T)*(1-np.cos(phi(r, v)))
        integrandi = lambda r, v: r*v*maxwell(v, T)*np.sin(phi(r, v))
        r = quad(lambda v: quad(lambda r: integrandr(r, v), 2e-8, 1.5e-6)[0], 1e2, 1e6)[0]
        i = quad(lambda v: quad(lambda r: integrandi(r, v), 2e-8, 1.5e-6)[0], 1e2, 1e6)[0]
        return r + 1j*i

    def __str__(self):
        result = f"States: {self.name}\n"
        result += f"Intensity: {self.intensity}\n"
        result += f"Unperturbed frequency: {self.coefficients[0]}\n"
        for i in range(1, len(self.coefficients)):
            result += f"C_{6*i} = {self.coefficients[i]}\n"
        return result
