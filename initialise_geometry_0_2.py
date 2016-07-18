import matplotlib.pyplot as plt
import numpy as np
from Param_cylinders_di_shell_0_2 import Shell_parameters

def Draw_cyl(R1_cyl, R2_cyl, gap, x_e, pos):

    R0, R1, R2, g2, gap_centre = Shell_parameters(R1_cyl, R2_cyl, gap)

    phi = np.arange(0, 2*np.pi, 0.01)

    #vertical
    if pos==1:

        x_e0 = max([R1_cyl, R2_cyl])

        #Trajectory in cylinder frame
        y_e = np.arange(-1000e-9, 1000e-9, 0.1e-9)
        x_e = np.ones(np.size(y_e)) * (x_e0+x_e)

        #Trajectory in annulus frame
        y_e_a = np.imag(g2/(x_e+1j*y_e) + 1j*R0)
        x_e_a = np.real(g2/(x_e+1j*y_e) + 1j*R0)

    #horizontal
    elif pos==2:

        #Trajectory in cylinder frame
        y_e = np.ones(np.size(x_e)) * \
                (gap_centre + gap/2 + 2*R1_cyl + x_e)
        x_e = np.arange(-500e-9, 500e-9, 0.1e-9)
        y_e = y_e * np.ones(np.size(x_e))

        #Trajectory in annulus frame
        y_e_a = np.imag(g2/(x_e+1j*y_e) + 1j*R0)
        x_e_a = np.real(g2/(x_e+1j*y_e) + 1j*R0)

    #in gap
    else:

        y_e = (gap_centre+x_e)
        x_e = np.arange(-500e-9, 500e-9, 0.1e-9)
        y_e = y_e * np.ones(np.size(x_e))

        #Trajectory in cylinder frame


        #Trajectory in annulus frame
        y_e_a = np.imag(g2/(x_e+1j*y_e) + 1j*R0)
        x_e_a = np.real(g2/(x_e+1j*y_e) + 1j*R0)




    #Plot the geometry in the annulus frame
    fig1 = plt.figure(1)
    plt.clf()
    ax1 = plt.gca()

    ax1.patch.set_facecolor('blue') #color area corresponding to cylinders
    ax1.patch.set_alpha(0.5)
    ax1.fill_between(R2*np.cos(phi), -R2*np.sin(phi), R2*np.sin(phi), color= 'white')
    ax1.fill_between(R1*np.cos(phi), R1*np.sin(phi), -R1*np.sin(phi), color= 'blue', alpha=0.5)

    plt.plot(x_e_a, y_e_a, 'r')     # plot the electron's trajectory
    plt.axis('equal')
    plt.ylim([-R2-2, R2+2])
    plt.figure(1).canvas.draw()

    # Plot the geometry in the cylinder frame
    fig2 = plt.figure(2)
    plt.clf()
    #Draw the cylinders and electron path
    circle1 = plt.Circle((0,gap_centre+gap/2+R1_cyl), R1_cyl,color='b', alpha=0.5)
    circle2 = plt.Circle((0, gap_centre-gap/2-R2_cyl), R2_cyl, color='b', alpha=0.5)
    fig2.gca().add_artist(circle1)
    fig2.gca().add_artist(circle2)
    plt.plot(x_e, y_e, 'r')             #Plots the electron's trajectory

    #Setting correct axis limits and ratios
    plt.axis('equal')
    xlim_n = -max((R1_cyl, R2_cyl)) - 3e-9

    if pos==2:
        plt.ylim([-2*R2_cyl+gap_centre-gap, y_e[0]+0.4e-9])
    else:
        plt.ylim([-2*R2_cyl+gap_centre-gap, 2*R1_cyl+gap+gap_centre])

    plt.xlim(xlim_n, -xlim_n)
    plt.figure(2).canvas.draw()




