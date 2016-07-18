import Tkinter as tk
from App_position_variables_0_2 import *
from initialise_geometry_0_2 import Draw_cyl
from Plot_power_absorbed_cylinders_0_2 import Absorption
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg

import numpy as np
import matplotlib.pyplot as plt

class Dimer:

    def __init__(self, parent, **kwargs):

        self.parent = parent
        #tk.Frame.__init__(self, parent, **kwargs)

    def draw_cyl(self):

        #Have to change label of electron position for gap case
        if self.orientation.get()==3:
            self.epos_var.set('Distance from gap centre in nm')

        #Getting the system parameters
        R1_cyl = float(self.R1_cyl_entry.get())*1e-9
        R2_cyl = float(self.R2_cyl_entry.get())*1e-9
        gap = float(self.gap_entry.get())*1e-9
        epos = float(self.epos_entry.get())*1e-9

        #Draw the system into the two empty canvases
        Draw_cyl(R1_cyl, R2_cyl, gap, epos, self.orientation.get())


    def plot_power(self):

        self.top_Q = tk.Toplevel(self.parent)
        self.top_Q.minsize(600,550)
        self.top_Q.configure(background='black')
        self.top_Q.wm_title('Power absorption')

        self.fig_3 = plt.figure(3, figsize=(5.5,4.5), dpi=100, facecolor='w')

        canvas_3 = FigureCanvasTkAgg(self.fig_3, master=self.top_Q)
        toolbar_3 = NavigationToolbar2TkAgg(canvas_3, self.parent)
        canvas_3.get_tk_widget().place(x=20,y=5)
        canvas_3.get_tk_widget()
        toolbar_3.place(x=20,y=480)

        R1_cyl = float(self.R1_cyl_entry.get())*1e-9
        R2_cyl = float(self.R2_cyl_entry.get())*1e-9
        gap = float(self.gap_entry.get())*1e-9
        epos = float(self.epos_entry.get())*1e-9
        vel = float(self.evel_entry.get())

        Absorption(R1_cyl, R2_cyl, gap, epos, vel, self.orientation.get())

    def draw_gui(self):

        self.subframe = tk.Frame(self.parent, width=300, height=400)
        self.subframe.place(x=c_x, y=c_y+30)
        self.subframe.config(background = 'black')

        self.orientation = tk.IntVar()
        self.orientation.set(1)

        ##############Radiobuttons##############################
        self.vert_rbutton = tk.Radiobutton(self.parent, text='vertical',\
                                      variable=self.orientation, value=1,\
                                      command=self.draw_cyl)
        self.vert_rbutton.place(x=c_x, width=RB_w,\
                                y=c_y+dy+c_h)

        self.gap_rbutton = tk.Radiobutton(self.parent, text='gap',\
                                     variable=self.orientation, value=3,\
                                     command=self.draw_cyl)
        self.gap_rbutton.place(x=c_x+dx+RB_w,\
                               y=c_y+dy+c_h,\
                               width=RB_w)

        self.hori_rbutton = tk.Radiobutton(self.parent, text='horizontal',\
                                           variable=self.orientation,\
                                           value=2, command= self.draw_cyl)
        self.hori_rbutton.place(x=c_x+2*dx+2*RB_w, \
                                y=c_y+dy + c_h,\
                                width=RB_w)


        #####################Draw some controls##################
        ##################FOR GEOMETRY############################
        #text and entry field cylinder radius R1
        self.R1_cyl_entry = tk.Entry(self.parent)
        self.R1_cyl_entry.place(x=c_x, y=c_y+2*(c_h+dy), width=e_w, height=c_h)
        self.R1_cyl_entry.delete(0,50)
        self.R1_cyl_entry.insert(0, '5')

        R1_cyl_label = tk.Label(self.parent, text='Radius of 1st cylinder in nm')
        R1_cyl_label.place(x=(e_w+dx+c_x), y=c_y+2*(c_h+dy), height=c_h, width = l_w)

        #text and entry field cylinder radius R2
        self.R2_cyl_entry = tk.Entry(self.parent)
        self.R2_cyl_entry.place(x=c_x, y=c_y+3*(c_h+dy), width=e_w, height=c_h)
        self.R2_cyl_entry.delete(0,50)
        self.R2_cyl_entry.insert(0, '5')

        R2_cyl_label = tk.Label(self.parent, text='Radius of 2nd cylinder in nm')
        R2_cyl_label.place(x=(e_w+dx+c_x), y=c_y+3*(c_h+dy), height=c_h, width = l_w)

        #text and entry field for gap between cylinder
        self.gap_entry = tk.Entry(self.parent)
        self.gap_entry.place(x=c_x, y=c_y+4*(c_h+dy), width=e_w, height=c_h)
        self.gap_entry.delete(0,50)
        self.gap_entry.insert(0, '0.4')

        gap_label = tk.Label(self.parent, text='Gap between cylinders in nm')
        gap_label.place(x=(e_w+dx+c_x), y=c_y+4*(c_h+dy), height=c_h, width = l_w)

        #text and entry field for electrons' postion cylinder
        self.epos_entry = tk.Entry(self.parent)
        self.epos_entry.place(x=c_x, y=c_y+5*(c_h+dy), width=e_w, height=c_h)
        self.epos_entry.delete(0,50)
        self.epos_entry.insert(0, '1')

        self.epos_var = tk.StringVar()
        self.epos_var.set('Distance of electron to cylinder in nm')
        epos_label = tk.Label(self.parent, textvariable = self.epos_var)
        epos_label.place(x=(e_w+dx+c_x), y=c_y+5*(c_h+dy), height=c_h, width = l_w)

        #draw_button
        draw_button = tk.Button(self.parent, text="Draw the geometry",\
                                    command = self.draw_cyl)
        draw_button.place(x=c_x,y=c_y+6*(c_h+dy), width=c_w)
        draw_button.configure(highlightbackground = 'black')

        #######Power absorption controls##############################
        #entry field for electron speed
        self.evel_entry = tk.Entry(self.parent)
        self.evel_entry.place(x=c_x, y=2*c_y+7*(c_h+dy), width=e_w, height=c_h)
        self.evel_entry.delete(0,50)
        self.evel_entry.insert(0, '0.05')

        evel_label = tk.Label(self.parent, text='Electron speed as a fraction of c')
        evel_label.place(x = (e_w+dx+c_x), y=2*c_y+7*(c_h+dy), height=c_h, width=l_w)


        Q_button = tk.Button(self.parent, text="Plot power absorption",\
                                    command = self.plot_power)
        Q_button.place(x=c_x,y=2*c_y+8*(c_h+dy), width=c_w)
        Q_button.configure(highlightbackground = 'black')








