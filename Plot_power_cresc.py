#Plot the potential of the periodic surface with a dipole source
import numpy as np
import matplotlib.pyplot as plt

from Calc_power_cresc import Pow_abs_rad, \
                                     Pow_abs_rad_r,\
                                     Pow_abs_rad_hori,\
                                     Pow_sca_rad, Pow_sca_r,\
                                     Pow_sca_hori

def Absorption(R1_cyl, R2_cyl, inv, epos, sca, vel, orientation):

    phi = np.arange(0,2*np.pi,0.001)

    z_1 = sca / (R1_cyl*np.exp(1j*phi) - inv)
    z_2 = sca / (R2_cyl*np.exp(1j*phi) - inv)

    #Geometrical and physical constants in the cylinder frame
    c = 3e8
    Conv = 1.602e-19/6.626e-34*2*np.pi #Conversion from eV to SI-units
    omega = np.arange(0.01, 8, 0.01)
    c_e = vel*c  #electron velocity

    plt.figure(4)
    plt.clf()
    plt.subplot(111)

    Q = np.zeros(np.size(omega))
    if orientation==1:
        x_e0 = min(np.real(z_2))
        x_e = x_e0 + np.sign(x_e0)*epos

        for m in range(0,np.size(omega)):
            Q[m] = Pow_abs_rad(inv, x_e, c_e, sca, R1_cyl, R2_cyl, omega[m])

        plt.plot(omega, Q/1.6e-19)

    elif orientation==3:
        x_e0 = max(np.real(z_2))
        x_e = x_e0 + epos

        for m in range(0,np.size(omega)):
            Q[m] = Pow_abs_rad_r(inv, x_e, c_e, sca, R1_cyl, R2_cyl, omega[m])

        plt.plot(omega, Q/1.6e-19)

    else:
        y_e0 = max(np.imag(z_1))
        x_e = y_e0 + epos

        for m in range(0,np.size(omega)):
            Q[m] = Pow_abs_rad_hori(inv, x_e, c_e, sca, R1_cyl, R2_cyl, omega[m])

        plt.plot(omega, Q/1.6e-19)


    plt.yscale('log')
    plt.xlabel('$\omega/eV$')
    plt.ylabel('Q/eV')
    plt.gcf().tight_layout()
    plt.figure(4).canvas.draw()

def Scattering(R1_cyl, R2_cyl, inv, epos, sca, vel, orientation):

    phi = np.arange(0,2*np.pi,0.001)

    z_1 = sca / (R1_cyl*np.exp(1j*phi) - inv)
    z_2 = sca / (R2_cyl*np.exp(1j*phi) - inv)

    #Geometrical and physical constants in the cylinder frame
    c = 3e8
    Conv = 1.602e-19/6.626e-34*2*np.pi #Conversion from eV to SI-units
    omega = np.arange(0.01, 8, 0.01)
    c_e = vel*c  #electron velocity

    plt.figure(5)
    plt.clf()
    plt.subplot(111)

    S = np.zeros(np.size(omega))
    if orientation==1:
        x_e0 = min(np.real(z_2))
        x_e = x_e0 + np.sign(x_e0)*epos

        for m in range(0,np.size(omega)):
            S[m] = Pow_sca_rad(inv, x_e, c_e, sca, R1_cyl, R2_cyl, omega[m])

        plt.plot(omega, S/1.6e-19)

    elif orientation==3:
        x_e0 = max(np.real(z_2))
        x_e = x_e0 + epos

        for m in range(0,np.size(omega)):
            S[m] = Pow_sca_r(inv, x_e, c_e, sca, R1_cyl, R2_cyl, omega[m])

        plt.plot(omega, S/1.6e-19)

    else:
        y_e0 = max(np.imag(z_1))
        x_e = y_e0 + epos

        for m in range(0,np.size(omega)):
            S[m] = Pow_sca_hori(inv, x_e, c_e, sca, R1_cyl, R2_cyl, omega[m])

        plt.plot(omega, S/1.6e-19)


    plt.yscale('log')
    plt.xlabel('$\omega/eV$')
    plt.ylabel('Scattering/eV')
    plt.gcf().tight_layout()
    plt.figure(5).canvas.draw()




