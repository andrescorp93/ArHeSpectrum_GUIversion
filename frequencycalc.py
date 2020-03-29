import numpy as np
import os
from calculationmethods import angtosm


cmtos1 = 2.998E10 * 2 * np.pi
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


def compound_states(statei, statef):
    if (statei[:2] == "1s") and (statef[:2] == "2p"):
        return statei + "_" + statef
    elif (statef[:2] == "1s") and (statei[:2] == "2p"):
        return statef + "_" + statei
    elif statef[:2] == statei[:2]:
        if int(statei[2:]) >= int(statef[2:]):
            return statei + "_" + statef
        else:
            return statef + "_" + statei


with open("intens.txt", "r") as intens_file:
    with open("energy_2.txt", "r") as energy_file:
        intens_txt = [line.strip().split("\t") for line in intens_file.readlines()]
        energy_txt = [line.strip().split("\t") for line in energy_file.readlines()]
        energies = {energy_txt[0][i]: np.array([float(energy_txt[k][i]) for k in range(1, len(energy_txt))])
                    for i in range(len(energy_txt[0]))}
        subkeys = [k for k in energies.keys()][1:]
        intens = {s[0]: [float(w) for w in s[1:]] for s in intens_txt[2:]}
        frequencies = {"R": energies["R"] * angtosm}
        for k in intens.keys():
            for i in range(len(intens[k])):
                if intens[k][i] != 0:
                    if find_state_by_substate(states, k) != find_state_by_substate(states, subkeys[i]):
                        nk = k + "_" + subkeys[i]
                        bs = compound_states(find_state_by_substate(states, k),
                                             find_state_by_substate(states, subkeys[i]))
                        frequencies[nk] = {"Intensity": intens[k][i],
                                           "Frequency Profile": np.abs(energies[k] - energies[subkeys[i]]) * cmtos1,
                                           "Base States": bs}

        frequencies_by_states = {"R": frequencies["R"]}
        for k, v in frequencies.items():
            if k != "R":
                if v["Base States"] not in frequencies_by_states.keys():
                    frequencies_by_states[v["Base States"]] = {
                        k: {"Intensity": v["Intensity"], "Frequency Profile": v["Frequency Profile"]}
                    }
                else:
                    frequencies_by_states[v["Base States"]][k] = {"Intensity": v["Intensity"],
                                                                  "Frequency Profile": v["Frequency Profile"]}

        try:
            os.mkdir("results")
        except OSError:
            pass
        finally:
            os.chdir("results")

        for k, v in frequencies_by_states.items():
            if k != "R":
                filename = k + ".txt"
                open(filename, "w").close()
                with open(filename, "w") as res_file:
                    start_string = ""
                    for sk in v.keys():
                        start_string += sk + "\t"
                    res_file.write(start_string.strip() + "\n")
                    intens_string = "Intensities\t"
                    for sk, sv in v.items():
                        if k != "R":
                            intens_string += str(sv["Intensity"]) + "\t"
                    res_file.write(intens_string.strip() + "\n")
                    for i in range(len(frequencies["R"])):
                        result = str(frequencies["R"][i])
                        for sk, sv in v.items():
                            if sk != "R":
                                result += "\t" + str(sv["Frequency Profile"][i])
                        res_file.write(result + "\n")
                    res_file.close()
