import numpy as np
import matplotlib.pyplot as plt
from calculationmethods import *
from spectrum import *

with open("results/coefficients/1s5_2p9.txt") as text_file:
    lines = [l.strip().split("\t") for l in text_file.readlines()]
    names = lines[0][1:]
    intensities = np.array([float(n) for n in lines[2][1:]])
    temperatures = np.array([float(l[0][:-6]) for l in lines[3::2]])
    broads = np.array([[float(s) for s in l[1:]] for l in lines[3::2]])
    shifts = np.array([[float(s) for s in l[1:]] for l in lines[4::2]])
    ns = np.arange(5e17, 3.5e18, 5e17)
    total_broads = np.zeros(len(temperatures))
    total_shifts = np.zeros(len(temperatures))
    for j in range(len(temperatures)):
        bs = np.zeros(len(ns))
        ws = np.zeros(len(ns))
        for k in range(len(ns)):
            n = ns[k]
            for i in range(len(names)):
                if i == 0:
                    spec = Spectrum(names[i], intensities[i], temperatures[j], n, broads[j, i], shifts[j, i])
                    spec.calculate_spectrum()
                else:
                    t = Spectrum(names[i], intensities[i], temperatures[j], n, broads[j, i], shifts[j, i])
                    t.calculate_spectrum()
                    spec += t
            ws[k], bs[k] = spec.fit_spectrum()
        # print(np.polyfit(ns, bs, 1)[0])
        # print(np.polyfit(ns, ws, 1)[0])
        total_broads[j] = np.polyfit(ns, bs, 1)[0]
        total_shifts[j] = np.polyfit(ns, ws, 1)[0]
    plt.plot(temperatures, total_broads, label='broads')
    plt.plot(temperatures, total_shifts, label='shifts')
    plt.legend()
    plt.show()
