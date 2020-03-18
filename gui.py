from tkinter import *
import os
from curve import Curve
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


def load_state():
    text.delete(1.0, END)
    text.insert(END, curves[statebox.curselection()[0]])


def plot_phase_shift():
    load_state()
    c = curves[statebox.curselection()[0]]
    plt.plot(c.grid, c.frequency, 'ro')
    plt.plot(c.grid, c.interpolated())
    plt.show()

def calc_cross_section():
    v = np.arange(1e4, 1e7, 5e3)
    c = curves[statebox.curselection()[0]]
    s_re = np.zeros(len(v))
    s_im = np.zeros(len(v))
    for i in range(len(v)):
        s = c.sigma_calc(v[i])
        s_re[i] = np.real(s)
        s_im[i] = np.imag(s)
        result = f"v = {v[i]} cm/s; s_shift = {s_im[i]} cm^2; s_broad = {s_re[i]} cm^2\n"
        text.insert(END, result)
    plt.plot(v, s_re)
    plt.plot(v, s_im)
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
calc_cross_button = Button(calc_frame, text="Calculate cross section", command=calc_cross_section)
calc_cross_button.grid(row=0, column=1)
for text_file in os.listdir("results"):
    filebox.insert(END, text_file)
    

root.mainloop()