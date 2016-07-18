import matplotlib.pyplot as plt
import numpy as np

def Draw_ell(semi_major, semi_minor, pos):

    axes_ratio = 1.0*semi_minor/semi_major

    u0 = np.arctanh(axes_ratio)
    scale = np.sqrt(semi_major**2 - semi_minor**2)

    R1_cyl = semi_major/scale + np.sqrt((semi_major/scale)**2 - 1)
    R2_cyl = semi_major/scale - np.sqrt((semi_major/scale)**2 - 1)
    phi = np.arange(0,2*np.pi,0.001)
    

    #Trajectory in cylinder frame
    y_e = np.linspace(-100*semi_major, 100*semi_major, 1e4)
    x_e = np.ones(np.size(y_e)) * (pos+semi_major)
    z_e = (x_e + 1j*y_e)/scale

    #Trajectory in annulus frame
    y_e_a_1 = np.imag(z_e + np.sqrt(z_e**2 - 1))
    x_e_a_1 = np.real(z_e + np.sqrt(z_e**2 - 1))
    y_e_a_2 = np.imag(z_e - np.sqrt(z_e**2 - 1))
    x_e_a_2 = np.real(z_e - np.sqrt(z_e**2 - 1))

    #Plot the geometry in the annulus frame
    fig1 = plt.figure(1)
    plt.clf()
    ax1 = plt.gca()

    ax1.patch.set_facecolor('white') #color area corresponding to cylinders
    ax1.fill_between(R1_cyl*np.cos(phi), R1_cyl*np.sin(phi),\
                    -R1_cyl*np.sin(phi), color= 'blue', alpha=0.5)
    ax1.fill_between(R2_cyl*np.cos(phi), -R2_cyl*np.sin(phi),\
                    R2_cyl*np.sin(phi), color= 'white')
    ax1.plot(R2_cyl*np.cos(phi), -R2_cyl*np.sin(phi), 'b')
    ax1.plot(R1_cyl*np.cos(phi), R1_cyl*np.sin(phi), 'b')

    plt.plot(x_e_a_1, y_e_a_1, 'r')
    plt.plot(x_e_a_2, y_e_a_2, 'r')     # plot the electron's trajectory
    plt.axis('equal')
    plt.ylim([-R2_cyl-1, max(R2_cyl+1, x_e_a_1[0])])
    plt.figure(1).canvas.draw()

    #Plot the geometry in ellipse frame
    fig2 = plt.figure(2)
    plt.clf()
    ax2 = plt.gca()

    z_1 = semi_major*np.cos(phi) + 1j*semi_minor*np.sin(phi)

    ax2.patch.set_facecolor('white') #color area of the ellipse
    ax2.fill_between(np.real(z_1), -np.imag(z_1),\
                    np.imag(z_1), color= 'blue', alpha=0.5)

    plt.plot(x_e, y_e, 'r')
    plt.axis('equal')
    
    plt.ylim([min(np.imag(z_1)-.5*semi_minor), max(np.imag(z_1)+.5*semi_minor)])
    plt.figure(2).canvas.draw()

    plt.figure(3)
    plt.plot(x_e_a_1, y_e_a_1, 'r')
    plt.plot(x_e_a_2, y_e_a_2, 'r')






