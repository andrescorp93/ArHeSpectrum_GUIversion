from calculation_methods import *
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

    def sigma_broad_calc(self):
        pass

    def sigma_shift_calc(self):
        pass

    def broad_coefficient_calc(self):
        pass

    def shift_coefficient_calc(self):
        pass
