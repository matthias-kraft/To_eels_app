import Tkinter as tk
from App_position_variables_0_2 import *
from initialise_crescent_0_2 import Draw_cres
from Plot_power_cresc import Absorption, Scattering
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
import matplotlib.pyplot as plt

class crescent:

    def __init__(self, parent, **kwargs):

        self.parent = parent

    def draw_cres(self):

        #Getting the system parameters
        R1_cyl = float(self.R1_cyl_entry.get())
        R2_cyl = float(self.R2_cyl_entry.get())
        inv = float(self.inv_entry.get())
        sca = float(self.scale_entry.get())*1e-9
        epos = float(self.epos_entry.get())*1e-9

        #Draw the system into the two empty canvases
        Draw_cres(R1_cyl, R2_cyl, inv, epos, sca, self.orientation.get())


    def plot_power(self):

        self.top_Q = tk.Toplevel(self.parent)
        self.top_Q.minsize(600,550)
        self.top_Q.configure(background='black')
        self.top_Q.wm_title('Power absorption')

        self.fig_4 = plt.figure(4, figsize=(5.5,4.5), dpi=100, facecolor='w')

        canvas_4 = FigureCanvasTkAgg(self.fig_4, master=self.top_Q)
        toolbar_4 = NavigationToolbar2TkAgg(canvas_4, self.parent)
        canvas_4.get_tk_widget().place(x=20,y=5)
        canvas_4.get_tk_widget()
        toolbar_4.place(x=20,y=480)

        #Getting the system parameters
        R1_cyl = float(self.R1_cyl_entry.get())
        R2_cyl = float(self.R2_cyl_entry.get())
        inv = float(self.inv_entry.get())
        sca = float(self.scale_entry.get())*1e-9
        epos = float(self.epos_entry.get())*1e-9
        vel = float(self.evel_entry.get())

        Absorption(R1_cyl, R2_cyl, inv, epos, sca, vel, self.orientation.get())

    def plot_scattering(self):

        self.top_S = tk.Toplevel(self.parent)
        self.top_S.minsize(600,550)
        self.top_S.configure(background='black')
        self.top_S.wm_title('Power absorption')

        self.fig_5 = plt.figure(5, figsize=(5.5,4.5), dpi=100, facecolor='w')

        canvas_5 = FigureCanvasTkAgg(self.fig_5, master=self.top_S)
        toolbar_5 = NavigationToolbar2TkAgg(canvas_5, self.parent)
        canvas_5.get_tk_widget().place(x=20,y=5)
        canvas_5.get_tk_widget()
        toolbar_5.place(x=20,y=480)

        #Getting the system parameters
        R1_cyl = float(self.R1_cyl_entry.get())
        R2_cyl = float(self.R2_cyl_entry.get())
        inv = float(self.inv_entry.get())
        sca = float(self.scale_entry.get())*1e-9
        epos = float(self.epos_entry.get())*1e-9
        vel = float(self.evel_entry.get())

        Scattering(R1_cyl, R2_cyl, inv, epos, sca, vel, self.orientation.get())


    def draw_gui(self):

        self.subframe = tk.Frame(self.parent, width=300, height=400)
        self.subframe.place(x=c_x, y=c_y+30)
        self.subframe.config(background = 'black')

        self.orientation = tk.IntVar()
        self.orientation.set(1)

        ##############Radiobuttons##############################
        self.left_rbutton = tk.Radiobutton(self.parent, text='left',\
                                      variable=self.orientation, value=1,\
                                      command=self.draw_cres)
        self.left_rbutton.place(x=c_x, width=RB_w,\
                                y=c_y+dy+c_h)

        self.gap_rbutton = tk.Radiobutton(self.parent, text='right',\
                                     variable=self.orientation, value=3,\
                                     command=self.draw_cres)
        self.gap_rbutton.place(x=c_x+dx+RB_w,\
                               y=c_y+dy+c_h,\
                               width=RB_w)

        self.right_rbutton = tk.Radiobutton(self.parent, text='top',\
                                           variable=self.orientation,\
                                           value=2, command= self.draw_cres)
        self.right_rbutton.place(x=c_x+2*dx+2*RB_w, \
                                y=c_y+dy + c_h,\
                                width=RB_w)


        #####################Draw some controls##################
        ##################FOR GEOMETRY############################
        #text and entry field cylinder radius R1
        self.R1_cyl_entry = tk.Entry(self.parent)
        self.R1_cyl_entry.place(x=c_x, y=c_y+2*(c_h+dy), width=e_w, height=c_h)
        self.R1_cyl_entry.delete(0,50)
        self.R1_cyl_entry.insert(0, '3')

        R1_cyl_label = tk.Label(self.parent, text='Radius of inner circle')
        R1_cyl_label.place(x=(e_w+dx+c_x), y=c_y+2*(c_h+dy), height=c_h, width = l_w)

        #text and entry field cylinder radius R2
        self.R2_cyl_entry = tk.Entry(self.parent)
        self.R2_cyl_entry.place(x=c_x, y=c_y+3*(c_h+dy), width=e_w, height=c_h)
        self.R2_cyl_entry.delete(0,50)
        self.R2_cyl_entry.insert(0, '5')

        R2_cyl_label = tk.Label(self.parent, text='Radius of outer circle')
        R2_cyl_label.place(x=(e_w+dx+c_x), y=c_y+3*(c_h+dy), height=c_h, width = l_w)

        #text and entry field for gap between cylinder
        self.inv_entry = tk.Entry(self.parent)
        self.inv_entry.place(x=c_x, y=c_y+4*(c_h+dy), width=e_w, height=c_h)
        self.inv_entry.delete(0,50)
        self.inv_entry.insert(0, '2')

        inv_label = tk.Label(self.parent, text='Inversion point')
        inv_label.place(x=(e_w+dx+c_x), y=c_y+4*(c_h+dy), height=c_h, width = l_w)

        #text and entry field for electrons' postion cylinder
        self.epos_entry = tk.Entry(self.parent)
        self.epos_entry.place(x=c_x, y=c_y+5*(c_h+dy), width=e_w, height=c_h)
        self.epos_entry.delete(0,50)
        self.epos_entry.insert(0, '1')

        self.epos_var = tk.StringVar()
        self.epos_var.set('Distance of electron to particle in nm')
        epos_label = tk.Label(self.parent, textvariable = self.epos_var)
        epos_label.place(x=(e_w+dx+c_x), y=c_y+5*(c_h+dy), height=c_h, width = l_w)

        #text and entry field for scale factor of trans postion cylinder
        self.scale_entry = tk.Entry(self.parent)
        self.scale_entry.place(x=c_x, y=c_y+6*(c_h+dy), width=e_w, height=c_h)
        self.scale_entry.delete(0,50)
        self.scale_entry.insert(0, '1')

        self.scale_var = tk.StringVar()
        self.scale_var.set('Scale factor divided by 10^(-9)')
        scale_label = tk.Label(self.parent, textvariable = self.scale_var)
        scale_label.place(x=(e_w+dx+c_x), y=c_y+6*(c_h+dy), height=c_h, width = l_w)


        #draw_button
        draw_button = tk.Button(self.parent, text="Draw the geometry",\
                                    command = self.draw_cres)
        draw_button.place(x=c_x,y=c_y+7*(c_h+dy), width=c_w)
        draw_button.configure(highlightbackground = 'black')

        #######Power absorption controls##############################
        #entry field for electron speed
        self.evel_entry = tk.Entry(self.parent)
        self.evel_entry.place(x=c_x, y=2*c_y+8*(c_h+dy), width=e_w, height=c_h)
        self.evel_entry.delete(0,50)
        self.evel_entry.insert(0, '0.05')

        evel_label = tk.Label(self.parent, text='Electron speed as a fraction of c')
        evel_label.place(x = (e_w+dx+c_x), y=2*c_y+8*(c_h+dy), height=c_h, width=l_w)


        Q_button = tk.Button(self.parent, text="Plot power absorption",\
                                    command = self.plot_power)
        Q_button.place(x=c_x,y=2*c_y+9*(c_h+dy), width=c_w)
        Q_button.configure(highlightbackground = 'black')

        S_button = tk.Button(self.parent, text="Plot photon scattering",\
                                    command = self.plot_scattering)
        S_button.place(x=c_x,y=2*c_y+10*(c_h+dy), width=c_w)
        S_button.configure(highlightbackground = 'black')


