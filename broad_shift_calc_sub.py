import matplotlib.pyplot as plt
import os
from scipy.integrate import nquad, simps, dblquad
from calculation_methods import *

T = np.arange(300, 850, 50)

os.chdir("results")
with open("1s5_2p10.txt", "r") as freq_file:
    lines = [l.strip().split("\t") for l in freq_file.readlines()]
    frequencies = {"R": np.array([float(lines[i][0]) for i in range(2, len(lines))])}
    for i in range(len(lines[0])):
        frequencies[lines[0][i]] = {"Intensity": float(lines[1][i + 1]),
                                    "Frequency Profile": np.array([float(lines[j][i + 1]) for j in range(2, len(lines))])}

    spectruminfo = {}
    for k in frequencies.keys():
        if k != "R":
            w0, c = interpolation_coefficients(frequencies["R"], frequencies[k]["Frequency Profile"])
            spectruminfo[k] = {"Zero Frequency": w0,
                               "Coefficients": c,
                               "<i|D|f><f|D|i>": frequencies[k]["Intensity"] * 1e-36}
    
    aij = 0
    for v in spectruminfo.values():
        aij += einstein_coefficient(v["<i|D|f><f|D|i>"], v["Zero Frequency"])

    print(aij / 3)
