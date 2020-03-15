import matplotlib.pyplot as plt
import numpy as np


def extract(infile, outfile):
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


extract("energies.out", "energy.txt")
