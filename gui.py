from tkinter import *
import os


def print_file():
    text.delete(1.0, END)
    text_file = open("results/" + lbox.get(lbox.curselection()), "r")
    text.insert(1.0, text_file.read())


root = Tk()
lbox = Listbox(height=20)
lbox.grid(row=0, column=0)
button = Button(text="Open", command=print_file)
button.grid(row=1, column=0)
scroll = Scrollbar(command=lbox.yview)
scroll.grid(row=0, column=1, sticky=N+S)
text = Text()
text.grid(row=0, column=2)
lbox.config(yscrollcommand=scroll.set)
for text_file in os.listdir("results"):
    lbox.insert(END, text_file)
    

root.mainloop()