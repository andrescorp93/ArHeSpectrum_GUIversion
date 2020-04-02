import numpy as np
import matplotlib.pyplot as plt

cmtos1 = 2.998E10 * 2 * np.pi
angtocm = 1E-8


def approx_coeffs(x, y, order=4, n=10):
    x_vec = [x ** (-order*i) for i in range(n+1)]
    m = np.array([[np.dot(x_vec[i], x_vec[j]) for j in range(n+1)] for i in range(n+1)])
    b = np.array([np.dot(y, x_vec[i]) for i in range(n+1)])
    return np.linalg.solve(m, b)


def potential(r, c, order=4):
    x = [r ** (-order*i) for i in range(len(c))]
    return np.dot(c, x)


with open("intens.txt", "r") as intens_file:
    with open("energy_2.txt", "r") as energy_file:
        intens_txt = [line.strip().split("\t") for line in intens_file.readlines()]
        energy_txt = [line.strip().split("\t") for line in energy_file.readlines()]
        energies = {energy_txt[0][i]: np.array([float(energy_txt[k][i]) for k in range(1, len(energy_txt))])
                    for i in range(len(energy_txt[0]))}
        for k in energies.keys():
            if k != "R":
                coefficients = approx_coeffs(energies["R"] * angtocm, energies[k])
                probe = np.array([potential(r * angtocm, coefficients) for r in energies["R"]])
                plt.plot(energies["R"] * angtocm, energies[k], "ro")
                plt.plot(energies["R"] * angtocm, probe)
                plt.show()
        
        