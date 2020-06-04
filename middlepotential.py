import numpy as np
import matplotlib.pyplot as plt

cmtoerg = 1.987E-16
daltontog = 1.661E-24
mu = daltontog * 40/11

states = {"1s5": ["1A1", "1A2", "1B1", "1B2", "2A2"], "1s4": ["2A1", "2B1", "2B2"],
          "1s3": ["3A2"], "1s2": ["3A1", "3B2", "3B1"], "2p10": ["4B2", "4B1", "4A2"],
          "2p9": ["5A2", "4A1", "5B2", "5B1", "6B2", "6B1", "6A2"],
          "2p8": ["5A1", "7B1", "7B2", "6A1", "7A2"], "2p7": ["8B1", "8B2", "8A2"],
          "2p6": ["7A1", "9B1", "9B2", "9A2", "8A1"], "2p5": ["9A1"],
          "2p4": ["10B1", "10B2", "10A2"], "2p3": ["10A1", "11A2", "11B2", "11B1", "11A1"],
          "2p2": ["12A2", "12B2", "12B1"], "2p1": ["12A1"]}


def find_state_by_substate(states, substate):
    for k in states.keys():
        try:
            states[k].index(substate)
            return k
        except ValueError:
            pass


with open("intens.txt", "r") as intens_file:
    with open("energy_2.txt", "r") as energy_file:
        intens_txt = [line.strip().split("\t") for line in intens_file.readlines()]
        energy_txt = [line.strip().split("\t") for line in energy_file.readlines()]
        energies = {energy_txt[0][i]: np.array([float(energy_txt[k][i]) for k in range(1, len(energy_txt))])
                    for i in range(len(energy_txt[0]))}
        middle_potentials = {"R": energies["R"]}
        for k, v in states.items():
            middle_potentials[k] = np.zeros(len(middle_potentials["R"]))
            n = len(v)
            for s in v:
                middle_potentials[k] += energies[s]
            middle_potentials[k] /= n
            middle_potentials[k] -= middle_potentials[k][-1]
            middle_potentials[k] *= cmtoerg / mu
        velocities = {"R": energies["R"]}
        for k, v in middle_potentials.items():
            if k != "R":
                velocities[k] = np.array([np.sqrt(v2) if v2 > 0 else 0 for v2 in v])
                for j in range(len(velocities[k]) - 2, 0, -1):
                    if velocities[k][j] < velocities[k][j + 1]:
                        velocities[k][j] = velocities[k][j + 1] 
                plt.plot(middle_potentials["R"], v)
        plt.show()
        for k, v in velocities.items():
            if k != "R":
                plt.plot(velocities["R"], v)
        plt.show()
        with open("velocities.txt", "w") as result:
            result.write("\t".join(velocities.keys()) + "\n")
            for i in range(len(velocities["R"])):
                result.write("\t".join([str(velocities[k][i]) for k in velocities.keys()]) + "\n")
