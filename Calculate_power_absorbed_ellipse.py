import numpy as np
from scipy.special import iv as BesselI

gamma = 0.032
n = np.arange(0,50,1)
kmax = max(n) + 20
ep0 = 8.85e-12
epD = 1.0
c = 3e8
Conv = 1.602e-19/6.626e-34*2*np.pi #Conversion from eV to SI-units

def Pow_abs(x_e, c_e, g, R0, R1, omega, epM):
    """Calculate the power absorbed in an annulus
    with a 'circling' electron as exciting source inside and
    an electron moving on a slightly curved trajectory outside (vertical)
    The trajectory of the electron derives from a straight vertical
    trajectory in the ellipse frame.
    Output: Resistive losses as a function of omega"""

    # epM = 1-64/(omega*(omega+1j*gamma))
    # omega = omega*Conv
    a_n_s = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf

    for k_n in range(0,np.size(n)):

        a_n_s[k_n] = np.exp(-omega/c_e*x_e)/omega*BesselI(k_n,omega*g/c_e)

    #Calculate expansion coefficients as in ELS_ellipse_annulus.pdf
    #This is for the cosine terms
    c_n_c_m = -2*a_n_s/((epM-1)**2 * R0**(2*n) - (epM+1)**2 * R1**(2*n))\
            * ((epM-1)*(R0*R1)**(2*n) + (epM+1)*R1**(2*n))

    c_n_c_p = -2*a_n_s/((epM-1)**2 * R0**(2*n) - (epM+1)**2 * R1**(2*n))\
            * ((epM-1) + (epM+1)*R1**(2*n))


    #This is for the sin terms

    c_n_s_m = -2*a_n_s/((epM-1)**2 * R0**(2*n) - (epM+1)**2 * R1**(2*n))\
            * ((epM-1)*(R0*R1)**(2*n) - (epM+1)*R1**(2*n))

    c_n_s_p = -2*a_n_s/((epM-1)**2 * R0**(2*n) - (epM+1)**2 * R1**(2*n))\
            * (-(epM-1) + (epM+1)*R1**(2*n))

    Term1 = sum(n**2/(2*n-1) * (pow(abs(c_n_s_p),2)+pow(abs(c_n_c_p),2) )\
            * (R1**(2*n)-R0**(2*n)))
    Term2 = sum(n**2/(2*n+1) * (pow(abs(c_n_s_m),2)+pow(abs(c_n_c_m),2) )\
            * (R1**(-2*n)-R0**(-2*n)))

    return np.pi * ep0/2 * omega * np.imag(epM) * (Term1-Term2)

def Pow_abs_r(x_e, c_e, g, R0, R1, omega, epM):
    """Calculate the power absorbed in an annulus
    with a 'circling' electron as exciting source inside and
    an electron moving on a slightly curved trajectory outside (vertical)
    The trajectory of the electron derives from a straight vertical
    trajectory in the ellipse frame.
    Output: Resistive losses as a function of omega"""

    #epM = 1-64/(omega*(omega+1j*gamma))
    #omega = omega*Conv
    k0 = omega/3e8
    a_n_s = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf

    for k_n in range(0,np.size(n)):

        a_n_s[k_n] = np.exp(-omega/c_e*x_e)/omega*BesselI(k_n,omega*g/c_e)

    #Calculate expansion coefficients as in ELS_ellipse_annulus.pdf
    #This is for the cosine terms
    c_n_c_m = -2*a_n_s/((epM-1)**2 * R0**(2*n) - (epM+1)**2 * R1**(2*n))\
            * ((epM-1)*(R0*R1)**(2*n) + (epM+1)*R1**(2*n))

    c_n_c_p = -2*a_n_s/((epM-1)**2 * R0**(2*n) - (epM+1)**2 * R1**(2*n))\
            * ((epM-1) + (epM+1)*R1**(2*n))


    #This is for the sin terms

    c_n_s_m = -2*a_n_s/((epM-1)**2 * R0**(2*n) - (epM+1)**2 * R1**(2*n))\
            * ((epM-1)*(R0*R1)**(2*n) - (epM+1)*R1**(2*n))

    c_n_s_p = -2*a_n_s/((epM-1)**2 * R0**(2*n) - (epM+1)**2 * R1**(2*n))\
            * (-(epM-1) + (epM+1)*R1**(2*n))

    #Radiative reaction
    #Denominator C
    #k0 = 0
    C0 = +1j*np.pi*g**2*k0**2/16.#/R0**2

    C = R0**2 * R1**(-2) * (-(epM-1)**2 * (R0/R1)**2 + (epM+1)**2)\
      + C0 * ((epM**2-1) * ((R0/R1)**2-1) * (R0**2 + R1**(-2)))\
      + 1*C0**2 * ((epM+1)**2*(R0/R1)**2 - (epM-1)**2)


    #Changed dipole terms
    c_n_c_m[1] = 2 * a_n_s[1] * R0**2 / C\
               *( (C0*(1-epM) + (1+epM)*R1**(-2))\
                - R1**(-2)*(C0*(1+epM) + R0**2*(1-epM)))

    c_n_c_p[1] = 2 * a_n_s[1] * R1**(-2) / C\
               *( ((epM-1)*R0**2*R1**(-2) - C0*(epM+1)*R0**2)\
               + ((epM+1)*R0**2 - (epM-1)*C0))

    c_n_s_m[1] = 2 * a_n_s[1] * R0**2 / C\
               *( -(C0*(1-epM) + (1+epM)*R1**(-2))\
                - R1**(-2)*(C0*(1+epM) + R0**2*(1-epM)))

    c_n_s_p[1] = 2 * a_n_s[1] * R1**(-2) / C\
               *( -((epM-1)*R0**2*R1**(-2) - C0*(epM+1)*R0**2)\
               + ((epM+1)*R0**2 - (epM-1)*C0))


    Term1 = sum(n**2/(2*n-1) * (pow(abs(c_n_s_p),2)+pow(abs(c_n_c_p),2) )\
            * (R1**(2*n)-R0**(2*n)))
    Term2 = sum(n**2/(2*n+1) * (pow(abs(c_n_s_m),2)+pow(abs(c_n_c_m),2) )\
            * (R1**(-2*n)-R0**(-2*n)))

    return np.pi * ep0/2 * omega * np.imag(epM) * (Term1-Term2)


def Pow_sca_r(x_e, c_e, g, R0, R1, omega, epM):
    """Calculate the power scattered by an annulus
    with a 'circling' electron as exciting source inside and
    an electron moving on a slightly curved trajectory outside (vertical)
    The trajectory of the electron derives from a straight vertical
    trajectory in the ellipse frame.
    Output: Resistive losses as a function of omega"""

    #epM = 1-64/(omega*(omega+1j*gamma))
    #omega = omega*Conv
    k0 = omega/c

    gamma_abs = 1j* np.pi**2 * ep0 * g**2/8.0 * k0**2
    k_n = 1

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf

    a_n_s = np.exp(-omega/c_e*x_e)/omega*BesselI(1,omega*g/c_e)

    #Radiative reaction
    #Denominator C
    C0 = +1j*np.pi*g**2*k0**2/16.#/R0**2

    C = R0**2 * R1**(-2) * (-(epM-1)**2 * (R0/R1)**2 + (epM+1)**2)\
      + C0 * ((epM**2-1) * ((R0/R1)**2-1) * (R0**2 + R1**(-2)))\
      + 1*C0**2 * ((epM+1)**2*(R0/R1)**2 - (epM-1)**2)


    #Calculate expansion coefficients as in ELS_ellipse_annulus.pdf
    #This is for the cosine terms
    b_c = -(a_n_s/C*( (C0*(epM+1)**2*(R0/R1)**2 - C0*(epM-1)**2\
                      - (epM**2-1) * ((R0/R1)**2-1)/R1**2)\
                      - 4*epM*(R0/R1)**2))\
          - 1*a_n_s

    #This is for the sin terms
    b_s = -(a_n_s/C*( -(C0*(epM+1)**2*(R0/R1)**2 - C0*(epM-1)**2\
                      - (epM**2-1) * ((R0/R1)**2-1)/R1**2)\
                      - 4*epM*(R0/R1)**2))\
          - 1*a_n_s


    return omega/2 * np.imag(gamma_abs * (abs(b_c)**2 + abs(b_s)**2))

def Pow_sca(x_e, c_e, g, R0, R1, omega, epM):
    """Calculate the power scattered by an annulus
    with a 'circling' electron as exciting source inside and
    an electron moving on a slightly curved trajectory outside (vertical)
    The trajectory of the electron derives from a straight vertical
    trajectory in the ellipse frame.
    Output: Resistive losses as a function of omega"""

    # epM = 1-64/(omega*(omega+1j*gamma))
    # omega = omega*Conv
    k0 = omega/3e8

    gamma_abs = 1j* np.pi**2 * ep0 * g**2/8 * k0**2
    k_n = 1

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in ELS_slab_crescent.pdf

    a_n_s = np.exp(-omega/c_e*x_e)/omega*BesselI(1,omega*g/c_e)

    #Calculate expansion coefficients as in ELS_ellipse_annulus.pdf
    #This is for the cosine terms
    b_c = (a_n_s /((epM-1)**2 * R0**(2) - (epM+1)**2 * R1**(2))\
          *( (epM**2-1) * (R1**(2)-R0**(2))\
             - 4*epM * R1**(2) * R0**(2) ) * R0**(-2)) - 1*a_n_s

    #This is for the sin terms
    b_s = (a_n_s/((epM-1)**2 * R0**(2) - (epM+1)**2 * R1**(2))\
          *( -(epM**2-1) * (R1**(2)-R0**(2))\
             - 4*epM * R1**(2) * R0**(2) ) * R0**(-2)) - 1*a_n_s


    return omega/2 * np.imag(gamma_abs * (abs(b_c)**2 + abs(b_s)**2))


if __name__ == "__main__":
    print 'Supposed to be called as a function, not main module'




