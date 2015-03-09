#!/usr/bin/env python

from Tkinter import *

root = Tk()

minutes = Label(root, text="Minutes:")
minutes.pack(side=LEFT)

scale = Scale(root, from_=1, to=45, orient=HORIZONTAL, length=300)
scale.pack()

button = Button(root, text="Hello", command=quit)
button.pack()

root.mainloop()