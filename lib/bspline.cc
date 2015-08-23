#include "cxxheaders.h"
#include "bspline.h"

/* Compute the B-spline function on a uniform mesh
 * Argument:
 *  x -- position of the singularity
 *  P -- number of B-spline points
 *  iBgn -- start point where the B-spline function value is non-zero
 *  w -- B-spline function values, it should be mapped to the physical
 *          axis as w[0] -> iBgn, w[1] -> iBgn+1, and so on
 * Note:
 *  -- The order of B-spline function is (P-1)
 *  -- The dimension of w must be at least P */
void computeBSpline(double x, int P, int &iBgn, double *w)
{
    // Alloc temp arrays
    double *u = new double [P];

    iBgn = (int)floor(x) - (P - 1);
    u[0] = iBgn - (x - P);
    for (int i = 1; i < P; i++) u[i] = u[i-1] + 1;

    w[0] = 1.0;
    for (int i = 1; i < P; ++i) w[i] = 0.0;

    for (int pp = 1; pp < P; ++pp) {
        double ipp = 1.0/pp;

        for (int j = pp; j >= 1; j--)
            w[j] = ipp*(u[j]*w[j] + (pp + 1 - u[j])*w[j-1]);
        w[0] = ipp*u[0]*w[0];
    } 

    // Dealloc temp arrays
    delete [] u;
}
