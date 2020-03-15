import numpy as np
import re

with open("DM.out", "r") as tdm_file:
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

    dx = np.array(dx)
    dy = np.array(dy)
    dz = np.array(dz)

    intens = np.absolute(dx) ** 2 + np.absolute(dy) ** 2 + np.absolute(dz) ** 2
    intens = intens[0]

    np.set_printoptions(precision=6)
    outfile = "TDM2.txt"
    symmetries = {0: "A1", 1: "B1", 2: "B2", 3: "A2"}
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
