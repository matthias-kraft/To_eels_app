#Plot the power absorbed by two nearly touching cylinders
#when a line electron flies past
import numpy as np
import matplotlib.pyplot as plt

from Calculate_power_absorbed_cylinders_0_2 \
import Pow_abs_vert, Pow_abs_gap, Pow_abs_hori

from Param_cylinders_di_shell_0_2 import Shell_parameters

def Absorption(R1_cyl, R2_cyl, gap, epos, vel, pos):

    #Geometrical and physical constants in the cylinder frame
    c = 3e8
    Conv = 1.602e-19/6.626e-34*2*np.pi #Conversion from eV to SI-units
    omega = np.arange(0.01, 8, 0.01)
    c_e = vel*c  #electron velocity

    ###Get the geometrical parameters in the di-shell frame
    R0, R1, R2, g2, gap_centre = Shell_parameters(R1_cyl, R2_cyl, gap)

    x_e_v = epos     #electron position vertical
    x_e_h = gap_centre + gap/2. + 2*R1_cyl + epos #electron position hori
    x_e_gap = gap_centre + epos

    Q = np.zeros(np.size(omega))

    #cross_sec = 2 * (R1_cyl+R2_cyl) + gap #Physical cross section

    plt.figure(3)
    plt.clf()
    plt.subplot(111)

    if pos==1: #vertical
        for m in range(0,np.size(omega)):
            Q[m] = Pow_abs_vert(R0, R1, R2, g2, x_e_v, c_e, omega[m])
        plt.plot(omega, Q, 'r')

    elif pos==2: #horizontal
        for m in range(0,np.size(omega)):
            Q[m] = Pow_abs_hori(R0, R1, R2, g2, x_e_h, c_e, omega[m])
        plt.plot(omega, Q, 'r')

    else: #gap
        for m in range(0,np.size(omega)):
            Q[m] = Pow_abs_gap(R0, R1, R2, g2, x_e_gap, c_e, omega[m])
        plt.plot(omega, Q, 'r')

    plt.yscale('log')
    plt.xlabel('$\omega/eV$')
    plt.ylabel('Q/eV')
    plt.gcf().tight_layout()
    plt.figure(3).canvas.draw()




