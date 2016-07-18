import numpy as np
from math import factorial

gamma = 0.032
n = np.arange(0,30,1)
ep0 = 8.85e-12
epD = 1.0
c = 3e8
Conv = 1.602e-19/6.626e-34*2*np.pi #Conversion from eV to SI-units

def Pow_abs(x0, x_e, c_e, g, R0, R1, omega):
    """Calculate the power absorbed in an annulus
    with a 'circling' electron as exciting source inside.
    The trajectory of the electron derives from a straight
    trajectory in the non-concentric annulus frame.
    Output: Resistive losses as a function of omega"""

    epM = 1-64/(omega*(omega+1j*gamma))
    omega = omega*Conv
    R0_x0 = R0/x0
    R1_x0 = R1/x0
    a_n_s = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf
    a_n_s[0] = np.exp(+omega/c_e*x_e) / omega

    for k_n in range(1,np.size(n)):
        k_s = n[1:k_n+1]
        sum_k = 0
        sum_k_in = 0
        for m in k_s:
            sum_k += factorial(k_n-1) * (-omega*g/c_e/x0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))
        a_n_s[k_n] = np.exp(+omega/c_e*x_e) / omega * sum_k
    #Calculate expansion coefficients as in ELS_slab_crescent.pdf
    b_n = a_n_s/\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (epD-epM) * (epD+epM) * (1.-(R1/R0)**(2*n))

    c_n = a_n_s/\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (-2.) * epD * (epD+epM) * R1_x0**(2*n)

    d_n = a_n_s/\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * 2. * epD * (epD-epM)

    Term1 = sum(n * pow(abs(c_n),2) * x0**(2*n) * (R0**(-2*n)-R1**(-2*n)))
    Term2 = sum(n * pow(abs(d_n),2) * x0**(-2*n) * (R1**(2*n)-R0**(2*n)))

    # dip_x = 2*sum(n * b_n * g / x0)
    # dip_y = -2j*sum(n * b_n * g / x0)

    return np.pi * ep0 * omega * np.imag(epM) * (Term1+Term2)
#    return dip_x, dip_y
#    return sum_k

def Pow_abs_rad_r(x0, x_e, c_e, g2, R0, R1, omega):
    """Calculate the power absorbed in an annulus
    with a 'circling' electron as exciting source inside.
    The trajectory of the electron derives from a straight
    trajectory in the non-concentric annulus frame, passing
    on the right of the annulus
    Output: Resistive losses as a function of omega"""

    epM = 1-64/(omega*(omega+1j*gamma))
    omega = omega*Conv
    k0 = omega/c
    R0_x0 = R0/x0
    R1_x0 = R1/x0
    a_n_s = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf
    a_n_s[0] = np.exp(+omega/c_e*x_e) / omega

    for k_n in range(1,np.size(n)):
        k_s = n[1:k_n+1]
        sum_k = 0
        sum_k_in = 0
        for m in k_s:
            sum_k += factorial(k_n-1) * (omega*g2/c_e/x0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))
        a_n_s[k_n] = np.exp(-omega/c_e*x_e) / omega * sum_k

    #Calculate scattering coefficients as in ELS_slab_crescent.pdf
    b_n_2 = 1./\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (epD-epM) * (epD+epM) * (1.-(R1/R0)**(2*n))

    #Calculate the radiative reaction term as in ELS_slab_radiative_reac
    a_n_r =  -1j*np.pi*(k0*g2/x0)**2/4  * np.sum(b_n_2*n*a_n_s)\
            / (1 + 1j*np.pi*(k0*g2/x0)**2/4 * np.sum(b_n_2*n))\
            * np.ones(np.size(n))

    # print a_n_r, a_n_s



    #Calculate expansion coefficients of induced field in metal
    c_n = (a_n_s + a_n_r)/\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (-2.) * epD * (epD+epM) * R1_x0**(2*n)

    d_n = (a_n_s + a_n_r)/\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * 2. * epD * (epD-epM)

    Term1 = sum(n * pow(abs(c_n),2) * x0**(2*n) * (R0**(-2*n)-R1**(-2*n)))
    Term2 = sum(n * pow(abs(d_n),2) * x0**(-2*n) * (R1**(2*n)-R0**(2*n)))


    return np.pi * ep0 * omega * np.imag(epM) * (Term1+Term2)


def Pow_abs_rad(x0, x_e, c_e, g2, R0, R1, omega):
    """Calculate the power absorbed in an annulus
    with a 'circling' electron as exciting source inside.
    The trajectory of the electron derives from a straight
    trajectory in the non-concentric annulus frame.
    Output: Resistive losses as a function of omega"""

    epM = 1-64/(omega*(omega+1j*gamma))
    omega = omega*Conv
    k0 = omega/c
    R0_x0 = R0/x0
    R1_x0 = R1/x0
    a_n_s = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf
    a_n_s[0] = np.exp(+omega/c_e*x_e) / omega

    for k_n in range(1,np.size(n)):
        k_s = n[1:k_n+1]
        sum_k = 0
        sum_k_in = 0
        for m in k_s:
            sum_k += factorial(k_n-1) * (-omega*g2/c_e/x0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))
        a_n_s[k_n] = np.exp(+omega/c_e*x_e) / omega * sum_k

    #Calculate scattering coefficients as in ELS_slab_crescent.pdf
    b_n_2 = 1./\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (epD-epM) * (epD+epM) * (1.-(R1/R0)**(2*n))

    #Calculate the radiative reaction term as in ELS_slab_radiative_reac
    a_n_r =  -1j*np.pi*(k0*g2/x0)**2/4  * np.sum(b_n_2*n*a_n_s)\
            / (1 + 1j*np.pi*(k0*g2/x0)**2/4 * np.sum(b_n_2*n))\
            * np.ones(np.size(n))

    # print a_n_r, a_n_s



    #Calculate expansion coefficients of induced field in metal
    c_n = (a_n_s + a_n_r)/\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (-2.) * epD * (epD+epM) * R1_x0**(2*n)

    d_n = (a_n_s + a_n_r)/\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * 2. * epD * (epD-epM)

    Term1 = sum(n * pow(abs(c_n),2) * x0**(2*n) * (R0**(-2*n)-R1**(-2*n)))
    Term2 = sum(n * pow(abs(d_n),2) * x0**(-2*n) * (R1**(2*n)-R0**(2*n)))


    return np.pi * ep0 * omega * np.imag(epM) * (Term1+Term2)

def Pow_sca_rad(x0, x_e, c_e, g2, R0, R1, omega):
    """Calculate the power scattered by the non-concentric
    annulus. This is equivalent to the power absorbed by a
    fictive absorver in the concentric annulus frame.
    Output: Power scattered by non-concentric annulus
            as a function of omega"""

    epM = 1-64/(omega*(omega+1j*gamma))
    omega = omega*Conv
    k0 = omega/c
    R0_x0 = R0/x0
    R1_x0 = R1/x0
    a_n_s = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf
    a_n_s[0] = np.exp(+omega/c_e*x_e) / omega

    for k_n in range(1,np.size(n)):
        k_s = n[1:k_n+1]
        sum_k = 0
        sum_k_in = 0
        for m in k_s:
            sum_k += factorial(k_n-1) * (-omega*g2/c_e/x0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))
        a_n_s[k_n] = np.exp(+omega/c_e*x_e) / omega * sum_k

    #Calculate scattering coefficients as in ELS_slab_crescent.pdf
    b_n_2 = 1./\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (epD-epM) * (epD+epM) * (1.-(R1/R0)**(2*n))

    #Calculate the radiative reaction term as in ELS_slab_radiative_reac
    a_n_r =  -1j*np.pi*(k0*g2/x0)**2/4  * np.sum(b_n_2*n*a_n_s)\
            / (1 + 1j*np.pi*(k0*g2/x0)**2/4 * np.sum(b_n_2*n))\
            * np.ones(np.size(n))

    #Calculate scattering coefficient including radiative losses

    b_n = b_n_2 * (a_n_r + a_n_s)

    E_sca_squ = 2. / x0**2 * abs(np.sum(b_n*n))**2

    return (np.pi*k0*g2)**2 * omega/4 * ep0 * E_sca_squ

def Pow_abs_rad_hori(x0, x_e, c_e, g2, R0, R1, omega):
    """Calculate the power absorbed in an annulus
    with a 'circling' electron as exciting source inside.
    The trajectory of the electron derives from a straight
    trajectory in the non-concentric annulus frame.
    Output: Resistive losses as a function of omega"""

    epM = 1-64/(omega*(omega+1j*gamma))
    omega = omega*Conv
    k0 = omega/c
    R0_x0 = R0/x0
    R1_x0 = R1/x0
    a_n_s = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf
    a_n_s[0] = np.exp(+omega/c_e*x_e) / omega

    for k_n in range(1,np.size(n)):
        k_s = n[1:k_n+1]
        sum_k = 0
        sum_k_in = 0
        for m in k_s:
            sum_k += factorial(k_n-1) * (1j*omega*g2/c_e/x0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))
        a_n_s[k_n] = np.exp(-omega/c_e*x_e) / omega * sum_k

    #Calculate scattering coefficients as in ELS_slab_crescent.pdf
    b_n_2 = 1./\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (epD-epM) * (epD+epM) * (1.-(R1/R0)**(2*n))

    #Calculate the radiative reaction term as in ELS_slab_radiative_reac
    a_n_r =  -1j*np.pi*(k0*g2/x0)**2/4  * np.sum(b_n_2*n*a_n_s)\
            / (1 + 1j*np.pi*(k0*g2/x0)**2/4 * np.sum(b_n_2*n))\
            * np.ones(np.size(n))

    # print a_n_r, a_n_s



    #Calculate expansion coefficients of induced field in metal
    c_n = (a_n_s + a_n_r)/\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (-2.) * epD * (epD+epM) * R1_x0**(2*n)

    d_n = (a_n_s + a_n_r)/\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * 2. * epD * (epD-epM)

    Term1 = sum(n * pow(abs(c_n),2) * x0**(2*n) * (R0**(-2*n)-R1**(-2*n)))
    Term2 = sum(n * pow(abs(d_n),2) * x0**(-2*n) * (R1**(2*n)-R0**(2*n)))


    return np.pi * ep0 * omega * np.imag(epM) * (Term1+Term2)

def Pow_sca_hori(x0, x_e, c_e, g2, R0, R1, omega):
    """Calculate the power scattered by the non-concentric
    annulus. This is equivalent to the power absorbed by a
    fictive absorver in the concentric annulus frame.
    Output: Power scattered by non-concentric annulus
            as a function of omega"""

    epM = 1-64/(omega*(omega+1j*gamma))
    omega = omega*Conv
    k0 = omega/c
    R0_x0 = R0/x0
    R1_x0 = R1/x0
    a_n_s = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf
    a_n_s[0] = np.exp(+omega/c_e*x_e) / omega

    for k_n in range(1,np.size(n)):
        k_s = n[1:k_n+1]
        sum_k = 0
        sum_k_in = 0
        for m in k_s:
            sum_k += factorial(k_n-1) * (1j*omega*g2/c_e/x0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))
        a_n_s[k_n] = np.exp(-omega/c_e*x_e) / omega * sum_k

    #Calculate scattering coefficients as in ELS_slab_crescent.pdf
    b_n_2 = 1./\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (epD-epM) * (epD+epM) * (1.-(R1/R0)**(2*n))

    #Calculate the radiative reaction term as in ELS_slab_radiative_reac
    a_n_r =  -1j*np.pi*(k0*g2/x0)**2/4  * np.sum(b_n_2*n*a_n_s)\
            / (1 + 1j*np.pi*(k0*g2/x0)**2/4 * np.sum(b_n_2*n))\
            * np.ones(np.size(n))

    #Calculate scattering coefficient including radiative losses

    b_n = b_n_2 * (a_n_r + a_n_s)

    E_sca_squ = 2. / x0**2 * abs(np.sum(b_n*n))**2
    return (np.pi*k0*g2)**2 * omega/4 * ep0 * E_sca_squ

def Pow_sca_r(x0, x_e, c_e, g2, R0, R1, omega):
    """Calculate the power scattered by the non-concentric
    annulus. This is equivalent to the power absorbed by a
    fictive absorver in the concentric annulus frame.
    Output: Power scattered by non-concentric annulus
            as a function of omega"""

    epM = 1-64/(omega*(omega+1j*gamma))
    omega = omega*Conv
    k0 = omega/c
    R0_x0 = R0/x0
    R1_x0 = R1/x0
    a_n_s = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf
    a_n_s[0] = np.exp(+omega/c_e*x_e) / omega

    for k_n in range(1,np.size(n)):
        k_s = n[1:k_n+1]
        sum_k = 0
        sum_k_in = 0
        for m in k_s:
            sum_k += factorial(k_n-1) * (omega*g2/c_e/x0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))
        a_n_s[k_n] = np.exp(-omega/c_e*x_e) / omega * sum_k

    #Calculate scattering coefficients as in ELS_slab_crescent.pdf
    b_n_2 = 1./\
            ((epD - epM)**2*R0_x0**(2*n) - (epD + epM)**2*R1_x0**(2*n))\
            * (epD-epM) * (epD+epM) * (1.-(R1/R0)**(2*n))

    #Calculate the radiative reaction term as in ELS_slab_radiative_reac
    a_n_r =  -1j*np.pi*(k0*g2/x0)**2/4  * np.sum(b_n_2*n*a_n_s)\
            / (1 + 1j*np.pi*(k0*g2/x0)**2/4 * np.sum(b_n_2*n))\
            * np.ones(np.size(n))


    #Calculate scattering coefficient including radiative losses

    b_n = b_n_2 * (a_n_r + a_n_s)

    E_sca_squ = 2. / x0**2 * abs(np.sum(b_n*n))**2

    return (np.pi*k0*g2)**2 * omega/4 * ep0 * E_sca_squ


if __name__ == "__main__":
    print 'Supposed to be called as a function, not main module'




