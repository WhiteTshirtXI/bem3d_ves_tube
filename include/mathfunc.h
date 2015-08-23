#ifndef MATHFUNC_H
#define MATHFUNC_H

#include "mkl.h"

// Some math constant
const double THRD = 1.0/3.0;

// In C++, the meaning of a%b for negative a is ambiguous.
// Here modulo(a,b) with b > 0 is always the NON-NEGATIVE remainder.
inline int modulo(int a, int b)
{
    assert(b > 0);
    int c = a%b;
    return (c >= 0)? c : c+b;
}


inline double modulo(double a, double b)
{
    assert(b > 0.0);
    return a - b*floor(a/b);
}


// square(x) = x*x
template<class T> 
T square(T x) 
{
    return x*x;
}


// Return a random number in [0, 1)
inline double rand01() 
{
    const int MAX = 32767;
    const double iMAX = 1.0/MAX;

    return iMAX*(rand()%MAX);
}


/* Determinant of a 2x2 matrix
 * Arguments:
 *  a -- the matrix */
inline double detMat2(const double (*a)[2])
{
    return a[0][0]*a[1][1] - a[0][1]*a[1][0];
}


/* Invert a 2x2 matrix
 * Arguments:
 *   a -- the original matrix
 *   b -- a^{-1} */
inline void invMat2(const double (*a)[2], double (*b)[2])
{
    double idetA = 1.0/detMat2(a);

    b[0][0] =  idetA*a[1][1];
    b[0][1] = -idetA*a[0][1];
    b[1][0] = -idetA*a[1][0];
    b[1][1] =  idetA*a[0][0];
}


/* Determinant of a 3x3 matrix
 * Arguments:
 *  a -- the matrix */
inline double detMat3(const double (*a)[3])
{
    return a[0][0]*a[1][1]*a[2][2] + a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1]
         - a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2] - a[0][2]*a[1][1]*a[2][0];
}


/* Invert a 3x3 matrix
 * Arguments:
 *   a -- the original matrix
 *   b -- a^{-1} */
inline void invMat3(const double (*a)[3], double (*b)[3])
{
    b[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
    b[0][1] = -(a[0][1]*a[2][2] - a[0][2]*a[2][1]);
    b[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];

    b[1][0] = -(a[1][0]*a[2][2] - a[1][2]*a[2][0]);
    b[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
    b[1][2] = -(a[0][0]*a[1][2] - a[0][2]*a[1][0]);

    b[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
    b[2][1] = -(a[0][0]*a[2][1] - a[0][1]*a[2][0]);
    b[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];

    double idetA = 1.0/detMat3(a);
    FOR_I3
    FOR_J3 {
        b[i][j] *= idetA;
    }
}


/* Evenly divide an interval [0, N) between processors
 * Arguments:
 *   N -- interval size
 *   nprocs -- number of processors 
 *   myrank -- my rank 
 *   ibgn -- beginning of local interval 
 *   iend -- end of local interval 
 *  Note:
 *   -- Each processor gets [ibgn, iend) */
inline void divideInterval(int N, int nprocs, int myrank, int &ibgn, int &iend)
{
    int m = N/nprocs;
    int r = modulo(N,nprocs);

    if (myrank < r) {
        ibgn = (m + 1)*myrank;
	iend = ibgn + (m + 1);
    } else {
        ibgn = m*myrank + r;
	iend = ibgn + m;
    }
}

#include "bspline.h"
#include "mblas.h"
#include "geom_oper.h"
#include "quadrature.h"

#endif
