from calculationmethods import *
from scipy.integrate import nquad, romberg, quad
from scipy.special import gamma

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

    def coefficient_calc(self, T, order=6):
        c = self.coefficients
        n = order
        p1 = np.power(2, (1/2)-(3/(n-1)))
        p2 = gamma(2-(1/(n-1))) * gamma((n-3)/(n-1))
        p3 = np.power(np.pi*mu/(R*T), -(1/2)+(1/(n-1)))
        p4 = np.power(-1j * gamma(n/2) / (gamma((n-1)/2) * c[1]), -2/(n-1))
        kappa = p1 * p2 * p3 * p4
        coeffs = np.zeros(len(c)-2)
        for j in range(2, len(c)):
            p1 = 1j * np.power(2, (5+n-2*j*n)/(1-n)) * np.power(np.pi, (3-j*n)/(2*(n-1)))
            p21 = gamma((n*j-3)/(n-1))
            p22 = gamma((n*j+1)/2)
            p23 = gamma((2*n*j+n-7)/(2*(n-1)))
            p24 = gamma((n*j-3)/(n-1)) * (n-1)
            p2 = p21*p22*p23 / p24
            p3 = np.power(mu/(R*T), (2+n-j*n)/(2*(n-1)))
            p4 = np.power(-1j * gamma(n/2) / (gamma((n-1)/2) * c[1]), (n*j-3)/(n-1))
            coeffs[j-2] = p1*p2*p3*p4
        return np.conj(kappa + np.dot(coeffs, c[2:]))

    def __str__(self):
        result = f"States: {self.name}\n"
        result += f"Intensity: {self.intensity}\n"
        result += f"Unperturbed frequency: {self.coefficients[0]}\n"
        for i in range(1, len(self.coefficients)):
            result += f"C_{6*i} = {self.coefficients[i]}\n"
        return result
