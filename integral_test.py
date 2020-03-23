from scipy.integrate import dblquad, nquad, simps, romberg, quad
from scipy.special import gamma
import numpy as np
from calculationmethods import phase_shift, maxwell
import matplotlib.pyplot as plt


p = np.array([7.74e-9, -8.54e-31, 2.34e-53, -3.09e-76, 1.97e-99, -4.90e-123])
t = 300

v = 5e4
n = 6
k = np.sqrt(np.pi)*gamma((n-1)/2)/gamma(n/2)
c = k * p[1] / v
eta = lambda r: r * np.sin(c * (r ** (-n+1)))
propose = (np.abs(c) ** (2/(n-1))) * gamma(-2/(n-1)) * np.sin(np.pi*(n-2)/(n-1)) / (n-1)
print(quad(eta, 0, np.Infinity))
print(propose)
# eta, 2e-8, np.Inf, lambda x: 0, lambda x: np.Infinity,

#print(nquad(eta, [[0,np.Infinity], [0,np.Infinity]],opts=[{'points': [[0,0],[0,np.Infinity],[np.Infinity,0]]},
#      {'epsabs': 1e-15},{'epsrel': 1e-8},{'limit': 1000}],full_output=True))

#print(romberg(lambda v: romberg(lambda r: eta(v, r), 2e-8, 1e-5), 1, 1e8))
