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
 * @file BSplineParamCurve.cpp
 * @author Alejandro Ungría Hirte
 * @version 1.0
 * @date 2020
 */

#include <assert.h>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <ios>
#include <cmath>

#include "BSpline.h"

#include "../build/FCMangle.h"

using namespace fitpackpp;

extern "C" {
	void parcur(int *iopt, int *ipar, int *dim, int *m, double *u, int *mx, double *x, double *w, double *ub, double *ue, int *k, double *s, int *nest, int *n, double *t, int *nc,double *c, double *fp, double *wrk, int *lwrk, int *iwrk, int *ier);
	void curev(int *idim, double *t, int *n, double *c, int *nc, int *k, double *u, int *m, double *x, int *mx, int *ier);
    void cualde(int *idim, double *t, int *n, double *c, int *nc, int *k1, double *u, double *d, int *nd, int *ier);
}

/**
 * @brief Constructor
 * @details Construct a B-Spline curve interpolation for the points specified by the list of abscissae x and ordinates y.
 * 
 * @param x Abscissae
 * @param y Ordinates
 * @param preferredDegree Preferred degree of the interpolating spline. 
 * The actual degree is chosen such as to be one less than the number of data points, but no higher than preferredDegree.
 * @param smoothing Smoothing factor.  Must be non-negative. Set to 0.0, i.e., no smoothing, by default.
 */
BSpline::BSpline(std::vector<double> &x, std::vector<double> &y, int preferredDegree, double smoothing)
{
    // std::cout << "We are here" << "\n";
	// Number of data points
	int m = (int) x.size();

    // std::cout << "The x array has a total of " << m << "members. \n";

    int idim = 2;

    int mx = (int) idim*m;

    

    std::vector<double> merged;

    auto v1 = x.begin ();
    auto v2 = y.begin ();

    while(v1 != x.end()){
        merged.push_back(*v1);
        merged.push_back(*v2);
        ++v1;
        ++v2;
    }

    // double* xv = &merged[0];

    // std::cout << merged.size() << "\n";

    
    

	// The actual degree of the spline must be less than m
	k = preferredDegree;
	if (k >= m) {
		k = m - 1;

		std::cerr << "WARNING:  Too few data points (" << m << ") to create B-Spline parametric curve of order " << preferredDegree << ". Reducing order to " << k << "." << std::endl;
	}

	// Configure curfit() parameters
	int iopt = 0;                       // Compute a smoothing spline
	int nest = (m + k + 1);               // Over-estimate the number of knots

	// Allocate weighting vector
	double *w = new double[m];
	std::fill(w, w + m, 1.0);

	// Allocate memory for knots and coefficients
	t = new double[nest];               // Knots
	std::fill(t, t + nest, 0.0);

    int nc = idim * nest;

	c = new double[nc];               // Coefficients
	std::fill(c, c + nc, 0.0);


	double fp = 0.0; // Weighted sum of squared residuals

	// Allocate working memory required by parcur
	int     lwrk = (m * (k + 1) + nest * (7 + idim + 3 * k));
	double *wrk  = new double[lwrk];
	std::fill(wrk, wrk + lwrk, 0.0);

	int    *iwrk = new int   [nest];
	std::fill(iwrk, iwrk + nest, 0);

    // Allocate working memory required by parcur
    u = new double[m];
    std::fill(u, u + m, 0.0);   
    // double u[m];

    int ipar = 0;
	int ier = 0;
    

    double ub = 0;
    double ue = 1;
    

    //void parcur(int *iopt, int *ipar, int *dim, int *m, double *u, int *mx, double *x, double *w, double *ub, double *ue, int *k, double *s, int *nest, int *n, double *t, double *c, double *fp, double *wrk, int *lwrk, int *iwrk, int *ier);
    //curfit(&iopt, &m, (double*) &x[0], (double*) &y[0], w, &x[0], &x[m - 1], &k, &smoothing, &nest, &n, t, c, &fp, wrk, &lwrk, iwrk, &ier);
	parcur(&iopt, &ipar, &idim, &m, u, &mx, (double*) &merged[0], w, &ub, &ue, &k, &smoothing, &nest, &n, t, &nc, c, &fp, wrk, &lwrk, iwrk, &ier);
	if (ier > 0) {
		if (ier >= 10) {
			std::stringstream s;
			s << "Error fitting B-Spline curve using parcur(): " << ier;
			throw std::runtime_error(s.str());
		} else {
			std::cerr << "WARNING:  Non-fatal error while fitting B-Spline curve using parcur(): " << ier << std::endl;
		}
	}
	std::vector<double> dist(x.size());
	double dx, dy;
	for(int i=0; i<x.size()-1; i++)
	{
		dx = x.at(i+1)-x.at(i);
		dy = y.at(i+1)-y.at(i);
		dist.at(i+1) = dist.at(i) + std::sqrt(dx*dx + dy*dy);
	}

	arc_length = dist.at(dist.size()-1);

    // std::cout << sizeof(n) << "\n";
	// De-allocate temporary memory
	delete[] w;
	delete[] wrk;
	delete[] iwrk;
    
}

BSpline::~BSpline(void)
{
	// Free memory
	delete[] t;
	delete[] c;
	delete[] u;
}

/**
 * @brief Knot X coordinates
 * @return Knot X coordinates
 */
std::vector<double> BSpline::knotX() 
{
	std::vector<double> knotX;
	knotX.assign(t, t + n);
	return knotX;
}

/**
 * @brief Coefficients
 * @return Coefficients
 */
std::vector<double> BSpline::coefs() 
{
	std::vector<double> coefs;
	coefs.assign(c, c + n);
	return coefs;
}	

/**
 * @brief Degree
 * @return Degree
 */
int BSpline::degree()
{
	return k;
}

std::vector<double> BSpline::curve_eval(std::vector<double> ts)
{

    int idim = 2;
	int m_pts = ts.size(); // Evaluate ts.size() points in total

	int ier = 0;
    int nc = coefs().size();

    // std::cout << "There are " << nc << " coefficients" << std::endl;
    int mx = idim*m_pts;

	std::vector<double> pts(mx);

    // std::cout  <<"There are " << mx << " entries."<<  "\n";
    // std::cout << "The degree is " << k << " \n";

	curev(&idim, t, &n, c, &nc, &k, (double *) &ts[0], &m_pts, (double *) &pts[0], &mx, &ier);
    
	if (ier > 0) {
		std::stringstream s;
		s << "Error evaluating B-Spline curve using curev() at point: " << ier;
		throw std::runtime_error(s.str());
	}

	return pts;
}


std::vector<double> BSpline::curv_der(double loc)
{
	
    int idim = 2;

    int ier = 0;
    int nc = coefs().size();
	
    int k1 = k + 1; //order

    int nd = k1 * idim;

	std::vector<double> d(nd);



	// splder also evaluates points in the exterior
	cualde(&idim, t, &n, c, &nc, &k1, &loc, (double *) &d[0],  &nd, &ier);
	if (ier > 0) {
		std::stringstream s;
		s << "Error evaluating B-Spline curve derivative using cualde() at point (" << loc <<  "): " << ier;
		throw std::runtime_error(s.str());
	}

	return d;

	
}

double BSpline::getCurvature(double s)
{
	std::vector<double> ds = curv_der(s);
	double dx = ds[2];
	double dy = ds[3];
	double dxdx = ds[4];
	double dydy = ds[5];

	double nominator = (dx*dydy - dy*dxdx);
	double denominator = (dx*dx+dy*dy);

	double kappa = nominator/(std::pow(denominator, 1.5));

	if(std::fabs(kappa) < 1e-5) kappa = 0;

	return kappa;


}

double BSpline::getAngle(double s)
{
	std::vector<double> derivatives = curv_der(s);
	double dx = derivatives[2];
	double dy = derivatives[3];

	double angle = std::atan2(dy,dx);

	return angle;

}