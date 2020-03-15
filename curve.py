from calculationmethods import *
class Curve:
    def __init__(self, grid, frequency, squared_dipole_moment):
        """Create a new state object"""
        self.frequency = frequency
        self.grid = grid
        self.interpolate()
        self.squared_dipole_moment = squared_dipole_moment
        self.intensity = einstein_coefficient(self.squared_dipole_moment, self.omega0)
    
    def interpolate(self):
        self.omega0, self.coefficients = interpolation_coefficients(self.grid, self.frequency)

    def eta(self, rho, s, x0, v):
        """
        Calculate eta for start position x0
        on the trajectory with length s
        on the distance rho from fixed pont
        with velocity v
        """
        return phase(rho, s, x0, v, self.coefficients)

    def sigma_calc(self):
        pass

    def coefficient_calc(self):
        pass

    def __str__(self):
        return super().__str__()
