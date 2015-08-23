#ifndef MBLAS_H
#define MBLAS_H

#include <algorithm>
#include "mkl.h"

/* double m_dmax 
 * double m_dmin
 * double m_dsum
 * double m_dnrm2
 * void m_dfill
 * void m_dclear
 * void m_dscal
 * double m_ddot
 * void m_dadd
 * void m_dsub
 * void m_daxpy
 * void m_dcopy
 * void m_dswap
 * void m_dgemv */

// max of a
inline double m_dmax(int N, const double *a)
{
    double r = -FLT_MAX;
    for (int i = 0; i < N; i++) r = (a[i] > r)? a[i] : r;
    return r;
}


// min of a
inline double m_dmin(int N, const double *a)
{
    double r = FLT_MAX;
    for (int i = 0; i < N; i++) r = (a[i] < r)? a[i] : r;
    return r;
}


// sum of a
inline double m_dsum(int N, const double *a)
{
    double sa = 0.0;
    for (int i = 0; i < N; i++) sa += a[i];
    return sa;
}

// L2 norm of a
inline double m_dnrm2(int N, const double *a)
{
    return cblas_dnrm2(N, a, 1);
}

// a[:] = s
inline void m_dfill(int N, double s, double *a)
{
    std::fill_n(a, N, s);
}


// a[:] = 0
inline void m_dclear(int N, double *a)
{
    std::fill_n(a, N, 0.0);
}


// a[:] = s * a[:]
inline void m_dscal(int N, double s, double *a)
{
    cblas_dscal(N, s, a, 1);
}


// dot product
inline double m_ddot(int N, const double *a, const double *b)
{
    return cblas_ddot(N, a, 1, b, 1);
}


// y = x
inline void m_dcopy(int N, const double *x, double *y)
{
    cblas_dcopy(N, x, 1, y, 1);
}


// y = y + x
inline void m_dadd(int N, const double *x, double *y)
{
    for (int i = 0; i < N; i++) y[i] += x[i];
}


// y = y - x 
inline void m_dsub(int N, const double *x, double *y)
{
    for (int i = 0; i < N; i++) y[i] -= x[i];
}


// y = y + alpha*x
inline void m_daxpy(int N, double alpha, const double *x, double *y)
{
    cblas_daxpy(N, alpha, x, 1, y, 1);
}


// swap x and y
inline void m_dswap(int N, double *x, double *y)
{
    cblas_dswap(N, x, 1, y, 1);
}


/* y = alpha * A * x + y
 * Arguments:
 *  m, n -- rows and columns
 *  A -- the matrix
 *  x -- 
 *  y -- */
inline void m_dgemv(int M, int N, double alpha, const double *A, const double *x, double *y)
{
    cblas_dgemv(CblasRowMajor, CblasNoTrans, M, N, alpha, A, N, x, 1, 1.0, y, 1);
}

#endif
