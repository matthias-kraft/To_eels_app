import numpy as np
from math import factorial

gamma = 0.32
n = np.arange(0,30,1)
ep0 = 8.85e-12
epD = 1.0
c = 3e8
Conv = 1.602e-19/6.626e-34*2*np.pi #Conversion from eV to SI-units

def Pow_abs_vert(R0, R1, R2, g2,  x_e, c_e, omega):
    """Calculate the power absorbed for two nearly
    touching cylinders if an electron passes them
    along their vertical axis.
    Output: Resistive losses as a function of omega
    Variable:
        R0 = inversion point
        R1 = inner radius of dielectric shell
        R2 = outer radius of dielectric shell
        g2  =   scaling factor g^2
        x_e = electrons position
        c_e = electron velocity
        omega = angular frequency in eV

    """

    epM = 1-64/(omega*(omega+1j*gamma))
    epD = 1.
    omega = omega*Conv
    beta = 2 * epD / (epM+epD)
    alpha = (epM-epD) / (epM+epD)

    an_p = np.zeros(np.size(n), dtype='complex128')
    an_n = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in Yu's document
    ###on EELS in nearly touching cylinders

    for k_n in range(1,np.size(n)):
        k_s = n[1:k_n+1]
        sum_k_n = 0
        sum_k_p = 0
        for m in k_s:
            sum_k_p += factorial(k_n-1) * (1j*omega*g2/c_e/R0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))
            sum_k_n += factorial(k_n-1) * (-1j*omega*g2/c_e/R0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))

        an_p[k_n] = (-1j)**k_n * np.exp(omega/c_e*(1j*g2/R0-x_e))\
                     / omega * sum_k_p
        an_n[k_n] = (1j)**k_n * np.exp(-omega/c_e*x_e)\
                     / omega * sum_k_n

    #Calculate expansion coefficients as in Yu's document

    fin_p = - beta * 1./(alpha**2-(R2/R1)**(2*n))\
            * (R2/R1)**(2*n) * an_p

    fin_n = beta * 1./(alpha**2-(R2/R1)**(2*n))\
            * alpha * (R0/R1)**(2*n) * an_n

    fon_p = beta * 1./(alpha**2-(R2/R1)**(2*n))\
            * alpha * (R2/R0)**(2*n) * an_p

    fon_n = - beta * 1./(alpha**2-(R2/R1)**(2*n))\
            * (R2/R1)**(2*n) * an_n

    #Calculate the power absorption

    Term1 = sum(n * (R1/R0)**(2*n) * (abs(fin_n)**2 + abs(fin_p)**2))
    Term2 = sum(n * (R0/R2)**(2*n) * (abs(fon_n)**2 + abs(fon_p)**2))

    return np.pi * ep0 * omega * np.imag(epM) * (Term1+Term2)

def Pow_abs_hori(R0, R1, R2, g2,  y_e, c_e, omega):
    """Calculate the power absorbed for two nearly
    touching cylinders if an electron passes them
    horizontally but outside the gap and not through
    the cylinders.
    Output: Resistive losses as a function of omega
    Variable:
        R0 = inversion point
        R1 = inner radius of dielectric shell
        R2 = outer radius of dielectric shell
        g2  =   scaling factor g^2
        y_e = electrons position vertical position
        c_e = electron velocity
        omega = angular frequency in eV

    """

    epM = 1-64/(omega*(omega+1j*gamma))
    epD = 1.
    omega = omega*Conv
    beta = 2 * epD / (epM+epD)
    alpha = (epM-epD) / (epM+epD)

    an_p = np.zeros(np.size(n), dtype='complex128')
    an_n = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in Yu's document
    ###on EELS in nearly touching cylinders

    for k_n in range(1,np.size(n)):
        k_s = n[1:k_n+1]
        sum_k_n = 0
        sum_k_p = 0
        for m in k_s:
            sum_k_p += factorial(k_n-1) * (-omega*g2/c_e/R0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))
            sum_k_n += factorial(k_n-1) * (+omega*g2/c_e/R0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))

        an_n[k_n] = (1j)**k_n * np.exp(omega/c_e*(g2/R0-y_e))\
                     / omega * sum_k_n
        an_p[k_n] = (-1j)**k_n * np.exp(-omega/c_e*y_e)\
                     / omega * sum_k_p


    #Calculate expansion coefficients as in Yu's document

    fin_n = - beta * 1./(alpha**2-(R2/R1)**(2*n))\
            * (R2/R1)**(2*n) * an_n

    fin_p = beta * 1./(alpha**2-(R2/R1)**(2*n))\
            * alpha * (R0/R1)**(2*n) * an_p

    fon_n = beta * 1./(alpha**2-(R2/R1)**(2*n))\
            * alpha * (R2/R0)**(2*n) * an_n

    fon_p = - beta * 1./(alpha**2-(R2/R1)**(2*n))\
            * (R2/R1)**(2*n) * an_p

    #Calculate the power absorption

    Term1 = sum(n * (R1/R0)**(2*n) * (abs(fin_n)**2 + abs(fin_p)**2))
    Term2 = sum(n * (R0/R2)**(2*n) * (abs(fon_n)**2 + abs(fon_p)**2))

    return np.pi * ep0 * omega * np.imag(epM) * (Term1+Term2)

def Pow_abs_gap(R0, R1, R2, g2,  y_e, c_e, omega):
    """Calculate the power absorbed for two nearly
    touching cylinders if an electron passes them
    horizontally exactly through the gap between the
    cylinders.
    Output: Resistive losses as a function of omega
    Variable:
        R0 = inversion point
        R1 = inner radius of dielectric shell
        R2 = outer radius of dielectric shell
        g2  =   scaling factor g^2
        y_e = electrons position vertical position
        c_e = electron velocity
        omega = angular frequency in eV

    """

    epM = 1-64/(omega*(omega+1j*gamma))
    epD = 1.
    omega = omega*Conv
    beta = 2 * epD / (epM+epD)
    alpha = (epM-epD) / (epM+epD)

    aS_p = np.zeros(np.size(n), dtype='complex128')
    aS_n = np.zeros(np.size(n), dtype='complex128')

    ###Lambda = 4*pi*eps0 in expression for source coefficients
    ###Calculate lambda according to formula in Yu's document
    ###on EELS in nearly touching cylinders

    for k_n in range(1,np.size(n)):
        k_s = n[1:k_n+1]
        sum_k_n = 0
        sum_k_p = 0
        for m in k_s:
            sum_k_p += factorial(k_n-1) * (-omega*g2/c_e/R0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))
            sum_k_n += factorial(k_n-1) * (-omega*g2/c_e/R0)**m\
                    / (factorial(m) * factorial(m-1) * \
                        factorial(k_n - m))

        aS_p[k_n] = (-1j)**k_n * np.exp(omega/c_e*(-g2/R0+y_e))\
                     / omega * sum_k_n
        aS_n[k_n] = (-1j)**k_n * np.exp(-omega/c_e*y_e)\
                     / omega * sum_k_p


    #Calculate expansion coefficients as in Yu's document

    bn_p = alpha * 1./(alpha**2-(R2/R1)**(2*n))\
            * ((R0/R1)**(2*n)*aS_n - alpha*aS_p)

    cn_p = alpha * 1./(alpha**2-(R2/R1)**(2*n))\
            * ((R2/R0)**(2*n)*aS_p - alpha*aS_n)

    fin_p = epD / epM * \
            (bn_p + aS_p - cn_p*(R0/R1)**(2*n))

    fon_p = epD / epM * \
            (cn_p + aS_n - bn_p*(R2/R0)**(2*n))

    #Calculate the power absorption

    Term1 = sum(n * (R1/R0)**(2*n) * abs(fin_p)**2)
    Term2 = sum(n * (R0/R2)**(2*n) * abs(fon_p)**2)

    return np.pi * ep0 * omega * np.imag(epM) * (Term1+Term2)


if __name__ == "__main__":
    print 'Error: Supposed to be called as a function, not main module.'




