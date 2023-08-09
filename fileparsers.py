import matplotlib.pyplot as plt
import numpy as np
import re
import os

c = 2.998E10 # speed of light
cmtos1 = c * 2 * np.pi
au_to_debye = 1 / 0.393456
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


def extract_energy(infile, outfile):
    with open(infile, "r") as en_file:
        m = [line.strip() for line in en_file.readlines() if line.strip() != '']
        en_string = [s for s in m if '0123456789'.find(s[0]) != -1]
        r_strings = [s[5:].replace("]", ", ") for s in m if s.find("R = [") != -1]
        r_string = ""
        for s in r_strings:
            r_string = r_string + s
        r_string = r_string[:-2]
        r = np.array([float(s) for s in r_string.split(", ")])
        a = []
        for s in en_string:
            t = s.split()
            a.append(float(t[2]))
        res = {}
        for i in range(len(a)):
            if (i % 48) not in res.keys():
                res[i % 48] = [a[i]]
            else:
                res[i % 48].append(a[i])
        for k in res.keys():
            res[k] = np.array(res[k])
        symmetries = {0: "A1", 1: "B1", 2: "B2", 3: "A2"}
        open(outfile, "w").close()
        with open(outfile, "w") as out:
            start_string = "R"
            for k in res.keys():
                start_string += "\t" + str(k % 12 + 1) + symmetries[k // 12]
            out.write(start_string + "\n")
            for i in range(len(r) - 1, -1, -1):
                s = str(r[i])
                for k in res.keys():
                    s += "\t" + str((res[k][i] - res[0][0]) * 219474.63 + 93143.7653)
                s += "\n"
                out.write(s)


def extract_dm(infile, outfile):
    with open(infile, "r") as tdm_file: # input from other file
        m = [line.strip() for line in tdm_file.readlines()]
        tdm_file.close()
        parsed_dx = []
        parsed_dy = []
        parsed_dz = []
        for j in range(len(m)):
            if re.search("DMX", m[j]) is not None:
                t1 = [s for s in m[j + 6: j + 884] if s != '']
                t = [s.split() for s in t1 if (re.search("Nr", s) is None) and (re.search("Total", s) is None)]
                for k in range(len(t)):
                    if k % 2 == 0:
                        t[k] = t[k][2:]
                parsed_dx.append(t)
            if re.search("DMY", m[j]) is not None:
                t1 = [s for s in m[j + 6: j + 884] if s != '']
                t = [s.split() for s in t1 if (re.search("Nr", s) is None) and (re.search("Total", s) is None)]
                for k in range(len(t)):
                    if k % 2 == 0:
                        t[k] = t[k][2:]
                parsed_dy.append(t)
            if re.search("DMZ", m[j]) is not None:
                t1 = [s for s in m[j + 6: j + 884] if s != '']
                t = [s.split() for s in t1 if (re.search("Nr", s) is None) and (re.search("Total", s) is None)]
                for k in range(len(t)):
                    if k % 2 == 0:
                        t[k] = t[k][2:]
                parsed_dz.append(t)
        dx = []
        dy = []
        dz = []
        for j in range(len(parsed_dx)):
            preout = []
            out = []
            for k in range(len(parsed_dx[j])):
                if k % 2 == 0:
                    r = np.array([float(s) for s in parsed_dx[j][k]])
                    i = np.array([float(s) for s in parsed_dx[j][k + 1]])
                    row = r + 1j * i
                    preout.append(row)
            for i in range(48):
                t1 = np.array(preout[i])
                t2 = np.array(preout[i + 48])
                t3 = np.array(preout[i + 96])
                t4 = np.array(preout[i + 144])
                t5 = np.array(preout[i + 192])
                t6 = np.array(preout[i + 240])
                out.append(np.hstack((t1, t2, t3, t4, t5, t6)))
            dx.append(np.array(out))
        for j in range(len(parsed_dy)):
            preout = []
            out = []
            for k in range(len(parsed_dy[j])):
                if k % 2 == 0:
                    r = np.array([float(s) for s in parsed_dy[j][k]])
                    i = np.array([float(s) for s in parsed_dy[j][k + 1]])
                    row = r + 1j * i
                    preout.append(row)
            for i in range(48):
                t1 = np.array(preout[i])
                t2 = np.array(preout[i + 48])
                t3 = np.array(preout[i + 96])
                t4 = np.array(preout[i + 144])
                t5 = np.array(preout[i + 192])
                t6 = np.array(preout[i + 240])
                out.append(np.hstack((t1, t2, t3, t4, t5, t6)))
            dy.append(np.array(out))
        for j in range(len(parsed_dz)):
            preout = []
            out = []
            for k in range(len(parsed_dz[j])):
                if k % 2 == 0:
                    r = np.array([float(s) for s in parsed_dz[j][k]])
                    i = np.array([float(s) for s in parsed_dz[j][k + 1]])
                    row = r + 1j * i
                    preout.append(row)
            for i in range(48):
                t1 = np.array(preout[i])
                t2 = np.array(preout[i + 48])
                t3 = np.array(preout[i + 96])
                t4 = np.array(preout[i + 144])
                t5 = np.array(preout[i + 192])
                t6 = np.array(preout[i + 240])
                out.append(np.hstack((t1, t2, t3, t4, t5, t6)))
            dz.append(np.array(out))
        dx = np.array(dx) * au_to_debye
        dy = np.array(dy) * au_to_debye
        dz = np.array(dz) * au_to_debye
        intens = np.absolute(dx) ** 2 + np.absolute(dy) ** 2 + np.absolute(dz) ** 2
        intens = intens[0]
        np.set_printoptions(precision=6)
        symmetries = ["A1", "B1", "B2", "A2"]
        open(outfile, "w").close()
        with open(outfile, "w") as out:
            start_string = "sqr(abs(dif))"
            for k in range(48):
                start_string += "\t" + str(k % 12 + 1) + symmetries[k // 12]
            out.write(start_string + "\n")
            for i in range(48):
                s = str(i % 12 + 1) + symmetries[i // 12]
                for k in range(i):
                    s += "\t" + str(intens[i][k])
                s += "\n"
                out.write(s)


def generate_freq_intense_profile(dmfilename, energyfilename, outdirname):
    with open(dmfilename, "r") as intens_file,  open(energyfilename, "r") as energy_file:
        intens_txt = [line.strip().split("\t") for line in intens_file.readlines()]
        energy_txt = [line.strip().split("\t") for line in energy_file.readlines()]
        energies = {energy_txt[0][i]: np.array([float(energy_txt[k][i]) for k in range(1, len(energy_txt))])
                    for i in range(len(energy_txt[0]))}
        limits = {}
        subkeys = [k for k in energies.keys()][1:]
        intens = {s[0]: [float(w) for w in s[1:]] for s in intens_txt[2:]}
        frequencies = {"R": energies["R"]}
        for k in intens.keys():
            for i in range(len(intens[k])):
                if intens[k][i] != 0:
                    if find_state_by_substate(states, k) != find_state_by_substate(states, subkeys[i]):
                        nk = k + "_" + subkeys[i]
                        bs = compound_states(find_state_by_substate(states, k),
                                             find_state_by_substate(states, subkeys[i]))
                        frequencies[nk] = {"Intensity": intens[k][i],
                                           "Frequency Profile": np.abs(energies[k] - energies[subkeys[i]])*cmtos1,
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
            os.mkdir(outdirname)
        except OSError:
            pass
        finally:
            os.chdir(outdirname)
        for k, v in frequencies_by_states.items():
            if k != "R":
                filename = k + ".txt"
                open(filename, "w").close()
                with open(filename, "w") as res_file:
                    start_string = "R\t"
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


def fullprepare(energyout="energies.out", dmout="DM.out", entxt="energy.txt", dm2txt="TDM2.txt", outdirname="results"):
    extract_energy(energyout, entxt)
    extract_dm(dmout, dm2txt)
    generate_freq_intense_profile(dm2txt, entxt, outdirname)
