# BSpline

This repository expands [Jorn Baayen's work on BSplines](https://github.com/jbaayen/fitpackpp) by including the possibility to smooth parametric curves. 

BSpline wraps the package of Fortran subroutines for smoothing splines by P. Dierckx, [FITPACK](http://www.netlib.org/dierckx/). This wrapper uses the double precision version of FITPACK distributed with [scipy](https://www.scipy.org). I

The BSplineCurve class wraps the 1D routines "curfit", "splev", and "splder".

The BSplineSurface class wraps the 2D routines "surfit", "bispev", and "parder".

The BSplineParamCurve class wraps the multidimensional routines "parcur", "curev" and "cualde".

# Installation

Pending.

# Usage

Pending.

# Testing 

Pendind.

# References

Paul Dierckx, Curve and Surface Fitting with Splines, Oxford University Press, 1993.

Jorn Baayen, C++ binding for FITPACK B-spline curve and surface fitting routines, [(GitHub repository)](https://github.com/jbaayen/fitpackpp)


