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
    for i in range(len(lines[0])):
        omega = np.array([float(lines[j][i + 1]) for j in range(2, len(lines))])
        curves.append(Curve(r, omega, float(lines[1][i + 1]) * 1e-36, names[i]))
    for curve in curves:
        statebox.insert(END, curve.name)

def calc_file():
    load_file()
    omega = [c.omega0 for c in curves]
    a = [c.intensity for c in curves]
    plt.plot(omega, a, 'ro')
    plt.show()


def load_state():
    text.delete(1.0, END)
    text.insert(END, curves[statebox.curselection()[0]])


def plot_phase_shift():
    load_state()
    c = curves[statebox.curselection()[0]]
    plt.plot(c.grid, c.frequency, 'ro')
    plt.plot(c.grid, c.interpolated())
    plt.show()


def calc_coefficients():
    load_state()
    c = curves[statebox.curselection()[0]]
    for t in np.arange(300, 1050, 50):
        r = c.coefficient_calc(t)
        print(r)
        text.insert(END, f"T = {t} K; ksi_b = {np.real(r)} cm^3/s; ksi_s = {np.imag(r)} cm^3/s\n")


def calc_subintegral():
    rho = np.arange(1.9e-8, 1e-7, 2e-9)
    T = 300
    c = curves[statebox.curselection()[0]]
    s = [c.sub_integrand(r, T) for r in rho]
    plt.plot(rho, np.real(s))
    plt.plot(rho, np.imag(s))
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
plot_button = Button(calc_frame, text="Plot", command=plot_phase_shift)
plot_button.grid(row=0, column=0)
calc_coeffs_button = Button(calc_frame, text="Calculate subintegral", command=calc_subintegral)
calc_coeffs_button.grid(row=0, column=1)
calc_coeffs_button = Button(calc_frame, text="Calculate coefficients", command=calc_coefficients)
calc_coeffs_button.grid(row=0, column=2)
for text_file in os.listdir("results"):
    filebox.insert(END, text_file)
    

root.mainloop()