import os
import numpy as np
from curve import Curve

for filename in os.listdir("results"):
    text_file = open("results/" + filename, "r")
    lines = [l.strip().split("\t") for l in text_file.readlines()]
    names = lines[0]
    curves = []
    r = np.array([float(lines[i][0]) for i in range(2, len(lines))])
    for i in range(len(lines[0])):
        omega = np.array([float(lines[j][i + 1]) for j in range(2, len(lines))])
        curves.append(Curve(r, omega, float(lines[1][i + 1]) * 1e-36, names[i]))
    einstein_total = sum([curve.intensity for curve in curves])
    print(f"Transition: {filename[:-4]} A_ki: {einstein_total}")
    
