# BSpline

BSpline wraps the package of Fortran subroutines for smoothing splines by P. Dierckx, FITPACK. This wrapper uses the double precision version of FITPACK distributed with scipy. It is a continuation of Jorn Baayen's Fitpackcpp wrapper: smoothing of parametric curves have been added.

The BSplineCurve class wraps the 1D routines "curfit", "splev", and "splder".

The BSplineSurface class wraps the 2D routines "surfit", "bispev", and "parder".

The BSplineParamCurve class wraps the n-D routines "parcur", "curev" and "cualde".

