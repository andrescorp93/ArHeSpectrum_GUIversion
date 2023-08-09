import numpy as np
from calculationmethods import *


class SpectrumError(Exception):
    def __init__(self, text):
        self.txt = text


class Spectrum:
    def __init__(self, name, intensity, temperature, n, kbroad=None, kshift=None, spec=None):
        self.name = name
        self.intensity = intensity
        self.temperature = temperature
        self.n = n
        self.kbroad = kbroad
        self.kshift = kshift
        self.spectrum = spec
        self.nu = np.arange(-1e10, 1.1e10, 1e7)

    def calculate_spectrum(self):
        delta = self.kshift * self.n
        width = self.kbroad * self.n
        self.spectrum = np.array([lorentz(f, self.intensity, delta, width) for f in self.nu])
        
    def fit_spectrum(self):
        ws = self.nu[np.argmax(self.spectrum)]
        m = np.max(self.spectrum)
        peakpart = []
        for i in range(len(self.nu)):
            if self.spectrum[i] > np.max(self.spectrum)/2:
                peakpart.append(self.nu[i])
        peakpart = np.array(peakpart)
        ds = (np.max(peakpart) - np.min(peakpart))
        a = np.pi * ds * m / 2
        popt, pcov = curve_fit(lorentz, self.nu, self.spectrum, p0=[a, ws, ds, 0])
        return popt[1], popt[2]

    def __add__(self, other):
        conditionstest = (self.n==other.n) and (self.temperature==other.temperature)
        try:
            if not conditionstest:
                raise SpectrumError("Not equal conditions!")
            newname = self.name + " with " + other.name
            newspectrum = self.spectrum + other.spectrum
            p = np.max(newspectrum)
            peakpart = []
            for i in range(len(self.nu)):
                if newspectrum[i] > p/2:
                    peakpart.append(self.nu[i])
            peakpart = np.array(peakpart)
            newwidth = np.max(peakpart) - np.min(peakpart)
            newintensity = np.pi*newwidth*p / 2
            return Spectrum(newname, newintensity, self.temperature, self.n, spec=newspectrum)
        except SpectrumError as sr:
            print(sr)
    
    def __iadd__(self, other):
        return self + other
