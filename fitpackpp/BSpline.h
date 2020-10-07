// Copyright (C) 2015 Deltares
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License version 2.1 as
// published by the Free Software Foundation.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

/**
 * @file BSpline.h
 * @author Alejandro Ungría Hirte
 * @version 1.0
 * @date 2020
 */

#ifndef BSPLINE_PARAM_CURVE_H
#define BSPLINE_PARAM_CURVE_H

#include <vector>
#include <string>

namespace fitpackpp
{

/**
 * @brief B-Spline curve
 * @details A wrapper for the curfit(), splev(), and splder() routines from FITPACK by  P. Dierckx:
 * http://www.netlib.org/fitpack/index.html
 * More specifically, we use the double-precision FITPACK version included with scipy:
 * https://github.com/scipy/scipy/tree/master/scipy/interpolate/fitpack
 */
class BSpline
{
public:
	BSpline(std::vector<double> &x, std::vector<double> &y, int preferredDegree=3, double smoothing=0.0);
	// BSpline(std::vector<double> &knotX, std::vector<double> &coefs, int degree);
	// BSpline(const std::string &filename);

	~BSpline();

	std::vector<double> knotX();
	std::vector<double> coefs();
	int degree();

    std::vector<double> curve_eval(std::vector<double> ts);
    std::vector<double> curv_der(double s);
	double getCurvature(double s);
	double getAngle(double s);

	double arc_length;
    // int *m;
    // int *mx;
    // int *idim;
    // int *nc;

private:
	int     k; // Spline degree
	int     n; // Number of knots
	double *t; // Knot coordinates
	double *c; // Spline coefficients
    double *u; // values u(i)
};

};

#endif