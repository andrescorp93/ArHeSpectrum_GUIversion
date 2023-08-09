from tkinter import *
import os
from curve import Curve
from calculationmethods import maxwell
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from spectrum import *
from fileparsers import *

curves = []


def ind(array, item):
    for idx, val in np.ndenumerate(array):
        if val == item:
            return idx


def tempcheck(value):
    newvalue = min(np.arange(300, 1050, 50), key=lambda x:abs(x-float(value)))
    temp_scale.set(newvalue)


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
    text.delete(1.0, END)
    for t in np.arange(300, 1050, 50):
        r = [curve.coefficient_calc(t) for curve in curves]
        text.insert(END, f'T={t} done!\n')
        text_file_out.write(str(t) + " Broad\t" +"\t".join([str(np.real(k)) for k in r]) + "\n")
        text_file_out.write(str(t) + " Shift\t" +"\t".join([str(np.imag(k)) for k in r]) + "\n")
    text_file_out.close()
    omega = [c.frequency[-2] for c in curves]
    a = [c.intensity for c in curves]
    fig.clear()
    spec_plot = fig.add_subplot()
    spec_plot.scatter(omega, a)
    canvas.draw()
    specfilebox.delete(0, END)
    for text_file in os.listdir("results/coefficients/"):
        if not os.path.isdir(os.path.join(os.path.abspath("results/coefficients/"),text_file)):
            specfilebox.insert(END, text_file)


def load_state():
    text.delete(1.0, END)
    text.insert(END, curves[statebox.curselection()[0]])


def plot_phase():
    load_state()
    c = curves[statebox.curselection()[0]]
    fig.clear()
    phase_plot = fig.add_subplot()
    phase_plot.plot(c.grid, c.phi)
    canvas.draw()


def calc_coefficients():
    load_state()
    c = curves[statebox.curselection()[0]]
    r = []
    text.insert(END, 'T, K\t broad, cm3/s\t shift, cm3/s\n')
    ts = np.arange(300, 1050, 50)
    for t in ts:
        r.append(c.coefficient_calc(t))
    for i in range(len(r)):
        text.insert(END, f'{ts[i]}\t{np.real(r[i]):.5e}\t{np.imag(r[i]):.5e}\n')
    fig.clear()
    coeff_plot = fig.add_subplot()
    coeff_plot.plot(np.arange(300, 1050, 50), np.real(r))
    coeff_plot.plot(np.arange(300, 1050, 50), np.imag(r))
    canvas.draw()
    
    
def calc_spectrum():
    text_file = open("results/coefficients/" + specfilebox.get(specfilebox.curselection()), "r")
    lines = [l.strip().split("\t") for l in text_file.readlines()]
    names = lines[0][1:]
    intensities = np.array([float(n) for n in lines[2][1:]])
    temperatures = np.array([float(l[0][:-6]) for l in lines[3::2]])
    broads = np.array([[float(s) for s in l[1:]] for l in lines[3::2]])
    shifts = np.array([[float(s) for s in l[1:]] for l in lines[4::2]])
    t = float(temp_scale.get())
    b = broads[ind(temperatures, t)]
    s = shifts[ind(temperatures, t)]
    n = float(conc_field.get()) if (conc_field.get() != '') else 1e18
    spec_data = Spectrum(names[0], intensities[0], t, n, b[0], s[0])
    spec_data.calculate_spectrum()
    for i in range(1, len(names)):
        tmp = Spectrum(names[i], intensities[i], t, n, b[i], s[i])
        tmp.calculate_spectrum()
        spec_data += tmp
    wfitted, bfitted = spec_data.fit_spectrum()
    text.insert(END, f'T = {t} K; n = {n} cm-3; dw = {wfitted:.5e} Hz; broad = {bfitted:.5e} Hz\n')
    fig.clear()
    spec_plot = fig.add_subplot()
    spec_plot.plot(spec_data.nu, spec_data.spectrum)
    # spec_plot.plot(spec_data.nu, spec_data.spectrum)
    canvas.draw()


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
calc_file_button = Button(file_frame, text="Calculate all transitions", command=calc_file)
calc_file_button.grid(row=2, column=0)

state_frame = Frame(root)
state_frame.grid(row=0, column=1)
statebox = Listbox(state_frame, height=30)
statebox.grid(row=0, column=0)
get_state_button = Button(state_frame, text="Load state", command=load_state)
get_state_button.grid(row=1, column=0)
scroll_state = Scrollbar(state_frame, command=statebox.yview)
scroll_state.grid(row=0, column=1, sticky=N+S)
statebox.config(yscrollcommand=scroll_state.set)

out_frame = Frame(root)
out_frame.grid(row=0, column=2)
text = Text(out_frame, height=30)
text.grid(row=0, column=0)
scroll_text = Scrollbar(out_frame, command=text.yview)
scroll_text.grid(row=0, column=1, sticky=N+S)

calc_frame = Frame(root)
calc_frame.grid(row=1, column=2)
plot_second_button = Button(calc_frame, text="Plot phase", command=plot_phase)
plot_second_button.grid(row=0, column=1)
calc_coeffs_button = Button(calc_frame, text="Calculate coefficients", command=calc_coefficients)
calc_coeffs_button.grid(row=0, column=2)

spec_frame = Frame(root)
spec_frame.grid(row=0, column=3)
spec_file_frame = Frame(spec_frame)
spec_file_frame.grid(row=0, column=1)
specfilebox = Listbox(spec_file_frame, height=20)
specfilebox.grid(row=0, column=0)
spec_param_frame = Frame(spec_frame)
spec_param_frame.grid(row=0, column=2)
temp_label = Label(spec_param_frame, text='Temperature')
temp_label.grid(row=0, column=1)
temp_scale = Scale(spec_param_frame, orient=HORIZONTAL, length=150, from_=300.0, to=1000, command=tempcheck)
temp_scale.grid(row=1, column=1)
conc_label = Label(spec_param_frame, text='Concentration')
conc_label.grid(row=2, column=1)
conc_field = Entry(spec_param_frame)
conc_field.grid(row=3, column=1)
spec_plot_frame = Frame(spec_param_frame)
spec_plot_frame.grid(row=4, column=1)
calc_coeffs_button = Button(spec_param_frame, text="Calculate spectrum", command=calc_spectrum)
calc_coeffs_button.grid(row=5, column=1)
fig = Figure()
canvas = FigureCanvasTkAgg(fig, master=spec_plot_frame)
canvas.draw()
canvas.get_tk_widget().pack()
toolbar = NavigationToolbar2Tk(canvas, spec_plot_frame)
toolbar.update()
canvas.get_tk_widget().pack()

for text_file in os.listdir('results'):
    if not os.path.isdir(os.path.join(os.path.abspath('results'),text_file)):
        filebox.insert(END, text_file)
try:
    for text_file in os.listdir("results/coefficients/"):
        if not os.path.isdir(os.path.join(os.path.abspath("results/coefficients/"),text_file)):
            specfilebox.insert(END, text_file)
except OSError:
    pass

root.mainloop()