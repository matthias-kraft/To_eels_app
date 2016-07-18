# Transformation optics for electron energy loss spectroscopy
This project provides a GUI interface to calculate electron energy loss and
cathodoluminescence spectra for a two-dimensional cylinder, a non-concentric
annulus and an ellipse. All calculations are based on a Transformation optics
approach to the problem that is presented in detail in:

1. [Matthias Kraft, Yu Luo and J.B. Pendry, Transformation Optics: A Time- and
Frequency-Domain Analysis of Electron-Energy Loss Spectroscopy, Nano Letters]
(http://pubs.acs.org/doi/abs/10.1021/acs.nanolett.6b02100)

2. [Yu Luo, Matthias Kraft and J.B. Pendry, Harnessing transformation optics for understanding electron energy loss and cathodoluminescence](https://arxiv.org/abs/1605.09319)

If this GUI app or the accompanying code helps you in your research, please
consider citing the two papers above.

I make no claims about the correctness of the code, as it is a modified version
of the one written for the two papers, so if you use it for your research please
double check that it works correctly.

# How to use the program
First, note that the calculations of the spectra are done in the quasi-static
limit. This means that the results are only valid for particles much smaller
than the wavelength of light at the plasmon resonance frequencies. As a rule
of thumb, our method is very accurate for particle sizes around 20nm and below
and provides good results for particle sizes up to 50-80nm. Also note that the
Transformation optics approach works in the non-relativistic limit, only. In
practice this means that the velocity of the exciting electron should be below
0.5c.

In all cases the permittivity of the material surrounding the nano-particles is
unity. The permittivity of the nano-particles is given by a Drude model,

<img src="http://mathurl.com/hsbaue5.png">.

Further note that all calculations are done in 2-D. That means all results are
for **line** electrons and in units per meter. For details read the papers.

In case you actually want to use the code provided here for you research, take
a look at the files:

  * Calculate_power_absorbed_cylinders_0_2.py
  * Plot_power_absorbed_cylinder_0_2.py
  * Calc_power_cresc.py
  * Plot_power_cresc.py
  * Calculate_power_absorbed_ellipse.py
  * Plot_power_ellipse.py

Unfortunately, the code itself is poorly documented at present, but feel free
to contact me with any questions you might have.

## Cylindrical dimer
Transformation optics works by relating a simple symmetric geometry to a less
symmetric one. It then exploits some properties of Maxwell's equation that allow
to transform the solution in the simple geometry, to the correct solution in the
complex geometry. For a cylindrical dimer the 2-D conformal transformation,

<img src="http://mathurl.com/zevulje.png">,
takes a ring-shaped void to a non-concentric annulus. In the program you can set
the following parameters:

  * The radius of the first cylinder.
  * The radius of the second cylinder.
  * The gap between the cylinders.
  * The distance of the electron to the dimer.
  * The orientation of the electron trajectory, i.e. whether it moves past the
  dimer **vertically**, **horizontally** or through the **gap** between the
  cylinders.
  * The electron's velocity

The program shows the original as well as the transformed system. The trajectory
of the electron is shown in red. Note that the electron has a constant velocity
in the cylindrical dimer frame but a non-constant one in the annulus frame. 

Press the **Plot power absorption** button to obtain the power absorbed by the
cylindrical dimer from the line electron. Note that this spectrum has been
calculated in the quasi-static limit.
electron.

## Non-concentric annulus
Much like the cylindrical dimer system, a non-concentric annulus can be
transformed to a concentric one by the conformal transformation,

<img src="http://mathurl.com/jvyxepp.png">.

In the program you can set the following parameters:

  * The radius of the inner circle of the concentric annulus.
  * The radius of the outer circle of the concentric annulus.
  * The inversion point x<sub>0</sub> of the transformation. This parameter
  **must** be **smaller** than the radius of the inner and outer circle.
  * The scale parameter g, which sets the size of the non-concentric annulus.
  * The distance of the electron to the non-concentric annulus.
  * The orientation of the electron trajectory, i.e. whether it moves past the
  non-concentric annulus **vertically to the left**, **vertically to the right**
  or **horizontally along the top**.
  * The electron's velocity

Press the **Plot power absorption**/**Plot power scattered** button to obtain
the power absorbed/scattered by the non-concentric annulus. Note that, while the
calculation has been carried out in the quasi-static limit, the effect
of radiation damping has been taken into account here.

## Ellipse
An ellipse too, can be transformed to a concentric annulus using the conformal
map,

<img src="http://mathurl.com/hs76jw6.png">.

In the program you can set the following parameters:

  * The semi-major axis of the ellipse.
  * The semi-minor axis of the ellipse.
  * The distance of the line electron to the ellipse
  * The electron's velocity

Press the **Plot power absorption**/**Plot power scattered** button to obtain
the power absorbed/scattered by the non-concentric annulus. Note that, while the
calculation has been carried out in the quasi-static limit, the effect
of radiation damping has been taken into account here.

# Installation

This program has been written in python 2.7.10. The GUI interface has been
written using TKinter, which is part of the standard python distribution.
However, the program also depends on Numpy, SciPy and Matplotlib. You need to
install these before you can run the program. If you are new to python the
easiest way to get all the necessary libraries is by installing the anaconda
package:

[https://www.continuum.io/downloads](https://www.continuum.io/downloads)
(please choose python 2.7.)

Once python is installed clone this git repository. Then open a terminal and
navigate (cd) inside the folder containing this program. In the folder, run 

    python TO_EELS_0_2.py

## Troubleshooting
If for some reason you are getting error messages when running the program as
described above, it may be best to install an exact copy of the packages that I
had installed into a [virtual environment]
(http://conda.pydata.org/docs/using/envs.html). The list of packages I used is given
in required_packages.txt. To setup a virtual environment with those packages
open a terminal and run
    
    conda create -n TO_EELS_app --file required_packages.txt

to activate the environment run

    source activate TO_EELS_app

with the virtual environment activated, navigate to the program's folder in your
terminal and again run

    python TO_EELS_0_2.py

This should work! To deactivate the virtual environment run

    source deactivate

If you still cannot run the program in a virtual environment and with the
packages and versions as in required_packages.txt, try google a solution. If you
are unsuccessful feel free to contact me.  

# MIT License

Copyright (c) 2016 Matthias Kraft, Yu Luo and J.B. Pendry

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.