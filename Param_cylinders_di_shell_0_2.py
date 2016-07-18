import numpy as np

def Shell_parameters(R1_cyl, R2_cyl, gap):
    """This calculates the inner and outer
    radius of the dielectric shell that one
    obtains from a shifted inverse of two
    nearly touching cylinders.
    Function returns the inversion point,
    inner Radius of shell,
    outer radius of shell and the scaling
    factor g^2, in this order.
    Input variables:
        R1_cyl = radius of first cylinder
        R2_cyl = radius of second cylinder
        gap = gap between the cylinders
    """

    ##Calculate position of inversion point in cylinder frame
    dp = (\
         -gap**2 - 2*gap*R1_cyl + \
         np.sqrt(gap) * np.sqrt(gap+2*R2_cyl) * \
         np.sqrt(gap**2 + 4*gap*R1_cyl + 4*R1_cyl**2 + \
              2*gap*R2_cyl + 4*R1_cyl*R2_cyl)\
         )\
         /(2*(gap + R1_cyl + R2_cyl))

    g = np.sqrt(dp) #np.sqrt of scaling factor

    d1 = g**2 / (2*R1_cyl+gap+dp) #inversion point to inner sphere
    d2 = g**2 / (2*R2_cyl-dp) #inversion point to outer sphere
    R1 = (g**2 / (gap+dp) - d1)/2. # inner radius of di-shell
    R2 = (g**2 / dp + d2)/2           #outer radius of di-shell
    R0 = R1 + d1         #Inversion point in di-shell frame

    #Calculate coordinate of the centre of the gap in the
    #cylinder frame
    gap_centre = - g**2*R0/(R2**2-R0**2) + R2_cyl + gap/2

    return R0, R1, R2, g**2, gap_centre

