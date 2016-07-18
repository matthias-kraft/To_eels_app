import matplotlib.pyplot as plt
import numpy as np

def Draw_cres(R1_cyl, R2_cyl, w0, x_e, g2, pos):

    phi = np.arange(0,2*np.pi,0.001)

    z_1 = g2 / (R1_cyl*np.exp(1j*phi) - w0)
    z_2 = g2 / (R2_cyl*np.exp(1j*phi) - w0)

    #left
    if pos==1:

        x_e0 = min(np.real(z_1))

        #Trajectory in cylinder frame
        y_e = np.linspace(-100*g2, 100*g2, 1e4)
        x_e = np.ones(np.size(y_e)) * (x_e0+np.sign(x_e0)*x_e)

        #Trajectory in annulus frame
        y_e_a = np.imag( g2 / (x_e+1j*y_e) + w0)
        x_e_a = np.real( g2 / (x_e+1j*y_e) + w0)

    #right
    elif pos==3:

        x_e0 = max(np.real(z_1))

        #Trajectory in cylinder frame
        y_e = np.linspace(-100*g2, 100*g2, 1e4)
        x_e = np.ones(np.size(y_e)) * (x_e0+x_e)

        #Trajectory in annulus frame
        y_e_a = np.imag( g2 / (x_e+1j*y_e) + w0)
        x_e_a = np.real( g2 / (x_e+1j*y_e) + w0)

    #top
    else:

        y_e0 = max(np.imag(z_1))

        #Trajectory in cylinder frame
        y_e = y_e0 + x_e
        x_e = np.arange(-500e-9, 500e-9, 0.1e-9)
        y_e = y_e * np.ones(np.size(x_e))

        #Trajectory in annulus frame
        y_e_a = np.imag( g2 / (x_e+1j*y_e) + w0)
        x_e_a = np.real( g2 / (x_e+1j*y_e) + w0)

    #Plot the geometry in the annulus frame
    fig1 = plt.figure(1)
    plt.clf()
    ax1 = plt.gca()

    ax1.patch.set_facecolor('white') #color area corresponding to cylinders
    ax1.fill_between(R2_cyl*np.cos(phi), -R2_cyl*np.sin(phi),\
                    R2_cyl*np.sin(phi), color= 'blue', alpha=0.5)
    ax1.fill_between(R1_cyl*np.cos(phi), R1_cyl*np.sin(phi),\
                    -R1_cyl*np.sin(phi), color= 'white')

    plt.plot(x_e_a, y_e_a, 'r')     # plot the electron's trajectory
    plt.axis('equal')
    plt.ylim([-R2_cyl-2, R2_cyl+2])
    plt.figure(1).canvas.draw()

    #Plot the geometry in non-concentric annulus frame
    fig2 = plt.figure(2)
    plt.clf()
    ax2 = plt.gca()

    ax2.patch.set_facecolor('white') #color area corresponding to cylinders
    ax2.fill_between(np.real(z_1), -np.imag(z_1),\
                    np.imag(z_1), color= 'blue', alpha=0.5)
    ax2.fill_between(np.real(z_2), np.imag(z_2),\
                    -np.imag(z_2), color= 'white')

    plt.plot(x_e, y_e, 'r')
    plt.axis('equal')

    if pos==2:
        plt.xlim([min(np.real(z_1)-1*g2), max(np.real(z_1)+1*g2)])
    else:
        plt.ylim([min(np.imag(z_1)-.5*g2), max(np.imag(z_1)+.5*g2)])
    plt.figure(2).canvas.draw()






