#Plot the power absorbed by two nearly touching cylinders
#when a line electron flies past
import numpy as np
import matplotlib.pyplot as plt

from Calculate_power_absorbed_ellipse \
import Pow_abs_r, Pow_sca_r


def Absorption(semi_major, semi_minor, epos, vel):

    #Geometrical and physical constants in the cylinder frame
    axes_ratio = 1.0*semi_minor/semi_major

    u0 = np.arctanh(axes_ratio)
    scale = np.sqrt(semi_major**2 - semi_minor**2)

    R1_cyl = semi_major/scale + np.sqrt((semi_major/scale)**2 - 1)
    R2_cyl = semi_major/scale - np.sqrt((semi_major/scale)**2 - 1)

    ##Use Drude model
    omega = np.arange(0.01,8,0.01)
    epM = 1 - 64/(omega*(omega + 1j*0.32))

    Conv = 1.602e-19/6.626e-34*2*np.pi #Conversion from eV to SI-units
    hbar = 1.0546e-34

    omega = omega*Conv

    #Geometrical and physical constants
    c = 3e8
    eV = 1.602e-19

    Q_r = np.zeros(np.size(omega))

    for m in range(0,np.size(omega)):
        Q_r[m] = Pow_abs_r(epos+semi_major, vel*c, scale, R2_cyl, R1_cyl,\
                            omega[m],epM[m])

    plt.figure(4)
    plt.clf()
    plt.subplot(111)    
    plt.plot(omega/Conv, Q_r/eV, 'r')
    plt.yscale('log')
    plt.xlabel('$\omega/eV$')
    plt.ylabel('Q/eV')
    plt.gcf().tight_layout()
    plt.figure(4).canvas.draw()

def Scattering(semi_major, semi_minor, epos, vel):

    #Geometrical and physical constants in the cylinder frame
    axes_ratio = 1.0*semi_minor/semi_major

    u0 = np.arctanh(axes_ratio)
    scale = np.sqrt(semi_major**2 - semi_minor**2)

    R1_cyl = semi_major/scale + np.sqrt((semi_major/scale)**2 - 1)
    R2_cyl = semi_major/scale - np.sqrt((semi_major/scale)**2 - 1)

    ##Use Drude model
    omega = np.arange(0.01,8,0.01)
    epM = 1 - 64/(omega*(omega + 1j*0.32))

    Conv = 1.602e-19/6.626e-34*2*np.pi #Conversion from eV to SI-units
    hbar = 1.0546e-34

    omega = omega*Conv

    #Geometrical and physical constants
    c = 3e8
    eV = 1.602e-19

    Sca_r = np.zeros(np.size(omega))

    for m in range(0,np.size(omega)):
        Sca_r[m] = Pow_sca_r(epos+semi_major, vel*c, scale, R2_cyl, R1_cyl,\
                             omega[m],epM[m])

    plt.figure(5)
    plt.clf()
    plt.subplot(111)    
    plt.plot(omega/Conv, Sca_r/eV, 'r')
    plt.yscale('log')
    plt.xlabel('$\omega/eV$')
    plt.ylabel('Pow$_{sca}$/eV')
    plt.gcf().tight_layout()
    plt.figure(5).canvas.draw()
    



