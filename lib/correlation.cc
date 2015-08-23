#include "cxxheaders.h"
#include "correlation.h"
#include "fftw3.h"
#include "mathfunc.h"

/* correlation
 * meanSqrDrift
 * correlation_direct
 * meanSqrDrift _direct*/


/* Calculate autocorrelation using FFT
 * Arguments:
 *   n -- number of samples 
 *   u -- samples
 *   ncorr -- number of correlation coefficients
 *   corr -- correlation coefficients
 * Note 
 *  -- corr(i) = sum_j=0^{n-1-i} u(j)*u(i+j) / (n-i) */
void correlation(int n, const double *u, int ncorr, double *corr)
{
    // Init
    for (int i = 0; i < ncorr; i++) corr[i] = 0.0;

    // Pad the input array with 0
    int n2 = 2;
    while (n2 < n+ncorr) n2*=2;

    double *uin = new double[n2];
    fftw_complex *uout = new fftw_complex[n2/2+1];

    fftw_plan fplan = fftw_plan_dft_r2c_1d(n2, uin, uout, FFTW_ESTIMATE);
    fftw_plan bplan = fftw_plan_dft_c2r_1d(n2, uout, uin, FFTW_ESTIMATE);

    for (int i = 0; i < n; i++) uin[i] = u[i];
    for (int i = n; i < n2; i++) uin[i] = 0.0;

    fftw_execute(fplan);
    for (int i = 0; i <= n2/2; i++) {
        double re = uout[i][0];
        double im = uout[i][1];
        uout[i][0] = re*re + im*im;
        uout[i][1] = 0.0;
    }
    fftw_execute(bplan);

    for (int i = 0; i < ncorr; i++) corr[i] = uin[i]/(n2*(n - i));

    // Dealloc temp arrays
    delete [] uin;
    delete [] uout;
    fftw_destroy_plan(fplan);
    fftw_destroy_plan(bplan);
} 


/* Calculate auto-correlation coefficient 
 * Arguments:
 *   n -- number of samples
 *   u -- samples
 *   ncorr -- number of correlation coefficients
 *   du2 -- mean square displacement
 * Note 
 *  -- du2(i) = \sum_j (u(j) - u(i+j))^2/(n-i)  0 <= j < n-i
 *            = \sum_j (u(j)^2 + u(i+j)^2 - 2*u(j)*u(i+j))/(n-i) */
void meanSqrDrift(int n, const double *u, int ncorr, double *du2)
{
    // Init
    for (int i = 0; i < ncorr; i++) du2[i] = 0.0;

    double *corr = new double[ncorr];
    correlation(n, u, ncorr, corr);

    du2[ncorr-1] = 0.0;
    for (int j = 0; j < n-(ncorr-1); j++) {
        du2[ncorr-1] += u[j]*u[j] + u[ncorr-1+j]*u[ncorr-1+j];
    }

    for (int i = ncorr-2; i >= 0; i--) {
        du2[i] = du2[i+1] + u[n-1-i]*u[n-1-i] + u[i]*u[i];
    }

    for (int i = 0; i < ncorr; i++) {
        du2[i] = du2[i]/(n-i) - 2*corr[i];
    }

    delete [] corr;
}


/* Directly calculate auto-correlation
 * Arguments:
 *   n -- number of samples
 *   u -- samples
 *   ncorr -- number of correlation coefficients
 *   corr -- correlation coefficients
 * Note 
 *  -- corr(i) = sum_j=0^{n-1-i} u(j)*u(i+j) / (n-i) */
void correlation_direct(int n, const double *u, int ncorr, double *corr)
{
    // Init
    for (int i = 0; i < ncorr; i++) corr[i] = 0.0;

    for (int i = 0; i < ncorr; i++) {
        int cnt = 0;
        for (int j = 0; j < n - i; j++) {
            corr[i] += u[j]*u[j+i];
            cnt++;
        }
        corr[i] /= cnt;
    }
}

/* Directly calculate the mean square drift
 * Arguments:
 *   n -- number of samples
 *   u -- samples
 *   ncorr -- number of correlation coefficiens
 *   du2 -- mean square displacement */
void meanSqrDrift_direct(int n, const double *u, int ncorr, double *du2)
{
    // Init
    for (int i = 0; i < ncorr; i++) du2[i] = 0.0;

    for (int i = 0; i < ncorr; i++) {
        int cnt = 0;
        for (int j = 0; j < n - i; j++) {
            du2[i] += square(u[j+i] - u[j]);
            cnt++;
        }
        du2[i] = du2[i]/cnt;
    }
}
