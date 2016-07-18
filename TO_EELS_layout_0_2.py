import Tkinter as tk
import numpy as np
import matplotlib.pyplot as plt
from App_position_variables_0_2 import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from initialise_geometry_0_2 import Draw_cyl
from Plot_power_absorbed_cylinders_0_2 import Absorption
import dimer_geom
import crescent_geom
import ellipse_geom

class TO_layout:

    def __init__(self, parent):

        self.parent = parent
        self.initialize()

    def draw_geom(self, parent):
        var = self.geometry.get()

        if var=="Cylindrical dimer":
            for widget in self.parent.winfo_children():
                widget.destroy()

            self.initialize()
            self.dimer.draw_gui()
            self.dimer.draw_cyl()


        elif var=="Non-concentric annulus":
            for widget in self.parent.winfo_children():
                widget.destroy()

            self.initialize()
            self.cresc.draw_gui()
            self.cresc.draw_cres()

        elif var=="Ellipse":
            for widget in self.parent.winfo_children():
                widget.destroy()

            self.initialize()
            self.ellipse.draw_gui()
            self.ellipse.draw_ellipse()



    def exit(self):
          self.parent.destroy()


    def initialize(self):

        self.frame = tk.Frame(self.parent, width=width, height=height)
        self.frame.configure(background='black')
        self.frame.pack()


        self.geometry = tk.StringVar(self.parent)
        self.geometry.set("Choose geometry") # default value

        select_geom = tk.OptionMenu(self.frame, self.geometry, "Cylindrical dimer",\
                                    "Non-concentric annulus", "Ellipse", command=self.draw_geom)
        select_geom.configure(background='black')
        select_geom.place(x=c_x, y=c_y, width=c_w)

        ####################Create exit button############################################
        exit = tk.Button(self.frame, text="Exit", command = self.exit)
        exit.place(x=c_x, y=c_y+fig_height*100, width=c_w)
        exit.configure(highlightbackground = 'black')

        #####################Create canvas to draw the geometry in #######################
        plt.close('all')
        self.fig_1 = plt.figure(1, figsize=(fig_width,fig_height), dpi=100, facecolor='w')

        canvas_1 = FigureCanvasTkAgg(self.fig_1, master=self.parent)
        canvas_1.get_tk_widget().place(x=fig1_x,y=fig1_y)
        canvas_1.get_tk_widget()

        self.fig_2 = plt.figure(2, figsize=(fig_width,fig_height), dpi=100, facecolor='w')

        canvas_2 = FigureCanvasTkAgg(self.fig_2, master=self.parent)
        canvas_2.get_tk_widget().place(x=fig2_x,y=fig2_y)
        canvas_2.get_tk_widget()

        ##################Create instances of all supported geometries####################
        self.dimer = dimer_geom.Dimer(self.parent)
        self.cresc = crescent_geom.crescent(self.parent)
        self.ellipse = ellipse_geom.ellipse(self.parent)










