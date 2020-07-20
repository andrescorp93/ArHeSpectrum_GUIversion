import os
import numpy as np
from curve import Curve
from calculationmethods import cmtos1

states = {"1s5": ["1A1", "1A2", "1B1", "1B2", "2A2"], "1s4": ["2A1", "2B1", "2B2"],
          "1s3": ["3A2"], "1s2": ["3A1", "3B2", "3B1"], "2p10": ["4B2", "4B1", "4A2"],
          "2p9": ["5A2", "4A1", "5B2", "5B1", "6B2", "6B1", "6A2"],
          "2p8": ["5A1", "7B1", "7B2", "6A1", "7A2"], "2p7": ["8B1", "8B2", "8A2"],
          "2p6": ["7A1", "9B1", "9B2", "9A2", "8A1"], "2p5": ["9A1"],
          "2p4": ["10B1", "10B2", "10A2"], "2p3": ["10A1", "11A2", "11B2", "11B1", "11A1"],
          "2p2": ["12A2", "12B2", "12B1"], "2p1": ["12A1"]}

for filename in os.listdir("results"):
    text_file = open("results/" + filename, "r")
    lines = [l.strip().split("\t") for l in text_file.readlines()]
    names = lines[0]
    curves = []
    r = np.array([float(lines[i][0]) for i in range(2, len(lines))])
    for i in range(len(lines[0]) - 1):
        omega = np.array([float(lines[j][i + 1]) for j in range(2, len(lines))])
        curves.append(Curve(r, omega, float(lines[1][i + 1]) * 1e-36, names[i]))
    g = len(states[filename[:-4].split("_")[1]])
    einstein_total = sum([curve.intensity for curve in curves])
    l = cmtos1 * 1e8 * len(curves) / sum([curve.frequency[-1] for curve in curves])
    print(f"Transition: {filename[:-4]} A_ki: {einstein_total/g} lambda: {l}")
    
