from tkinter import *
import os
from curve import Curve
from calculationmethods import maxwell
import numpy as np
import matplotlib.pyplot as plt

curves = []

def load_file():
    text.delete(1.0, END)
    statebox.delete(0, END)
    curves.clear()
    text_file = open("results/" + filebox.get(filebox.curselection()), "r")
    lines = [l.strip().split("\t") for l in text_file.readlines()]
    names = lines[0]
    r = np.array([float(lines[i][0]) for i in range(2, len(lines))])
    for i in range(len(lines[0])-1):
        omega = np.array([float(lines[j][i + 1]) for j in range(2, len(lines))])
        curves.append(Curve(r, omega, float(lines[1][i + 1]) * 1e-36, names[i + 1]))
    for curve in curves:
        statebox.insert(END, curve.name)

def calc_file():
    load_file()
    try:
        os.mkdir("results/coefficients/")
    except OSError:
        pass
    text_file_out = open("results/coefficients/" + filebox.get(filebox.curselection()), "w")
    text_file_out.write("States\t" +"\t".join([curve.name for curve in curves]) + "\n")
    text_file_out.write("Frequency\t" +"\t".join([str(curve.frequency[-1]) for curve in curves]) + "\n")
    text_file_out.write("Intensity\t" +"\t".join([str(curve.intensity) for curve in curves]) + "\n")
    for t in np.arange(300, 1050, 50):
        r = [curve.coefficient_calc(t) for curve in curves]
        text_file_out.write(str(t) + " Broad\t" +"\t".join([str(np.real(k)) for k in r]) + "\n")
        text_file_out.write(str(t) + " Shift\t" +"\t".join([str(np.imag(k)) for k in r]) + "\n")
    text_file_out.close()
    omega = [c.frequency[-2] for c in curves]
    a = [c.intensity for c in curves]
    plt.plot(omega, a, 'ro')
    plt.show()


def load_state():
    text.delete(1.0, END)
    text.insert(END, curves[statebox.curselection()[0]])


def plot_phase():
    load_state()
    c = curves[statebox.curselection()[0]]
    plt.plot(c.r, c.r*np.sin(c.phi/100000))
    plt.show()


def calc_coefficients():
    load_state()
    c = curves[statebox.curselection()[0]]
    r = []
    for t in np.arange(300, 1050, 50):
        r.append(c.coefficient_calc(t))
    for x in r:
        print(x)
    plt.plot(np.arange(300, 1050, 50), np.real(r))
    plt.plot(np.arange(300, 1050, 50), np.imag(r))
    plt.show()


root = Tk()
file_frame = Frame(root)
file_frame.grid(row=0, column=0)
filebox = Listbox(file_frame, height=20)
filebox.grid(row=0, column=0)
get_file_button = Button(file_frame, text="Open", command=load_file)
get_file_button.grid(row=1, column=0)
scroll_file = Scrollbar(file_frame, command=filebox.yview)
scroll_file.grid(row=0, column=1, sticky=N+S)
filebox.config(yscrollcommand=scroll_file.set)
calc_file_button = Button(file_frame, text="Calculate spectrum", command=calc_file)
calc_file_button.grid(row=2, column=0)

state_frame = Frame(root)
state_frame.grid(row=0, column=1)
statebox = Listbox(state_frame, height=20)
statebox.grid(row=0, column=0)
get_state_button = Button(state_frame, text="Load state", command=load_state)
get_state_button.grid(row=1, column=0)
scroll_state = Scrollbar(state_frame, command=statebox.yview)
scroll_state.grid(row=0, column=1, sticky=N+S)
statebox.config(yscrollcommand=scroll_state.set)

out_frame = Frame(root)
out_frame.grid(row=0, column=2)
text = Text(out_frame)
text.grid(row=0, column=0)
scroll_text = Scrollbar(out_frame, command=text.yview)
scroll_text.grid(row=0, column=1, sticky=N+S)

calc_frame = Frame(root)
calc_frame.grid(row=1, column=2)
plot_second_button = Button(calc_frame, text="Plot phase", command=plot_phase)
plot_second_button.grid(row=0, column=1)
calc_coeffs_button = Button(calc_frame, text="Calculate sigmas", command=calc_coefficients)
calc_coeffs_button.grid(row=0, column=2)
for text_file in os.listdir("results"):
    filebox.insert(END, text_file)
    

root.mainloop()