#ifndef HISTOGRAM_H
#define HISTOGRAM_H

/* Calculate the histogram
 * Arguments:
 *   u -- the sample
 *   n -- number of samples
 *   umin, umax -- the upper and lower bound
 *   nrho -- the number of intervals in the histogram
 *   rho -- the probability density
 * Note:
 *   -- The normalization is such that the integral of rho is 1 */
void histogram(const int n, const double *u, 
		const double umin, const double umax, 
		const int nrho, double *rho)
{
    const double du = (umax - umin)/nrho;
    const double idu = 1.0/du;

    // Init
    for (int i = 0; i < nrho; i++) rho[i] = 0.0;

    // Count number of samples in each interval
    for (int p = 0; p < n; p++) {
        int i = (int)floor( (u[p] - umin)*idu );

	if (i >= 0 && i < nrho) rho[i]++;
    }

    // Normalization
    double rho_sum = 0.0;
    for (int i = 0; i < nrho; i++) rho_sum += rho[i];

    double scal = 1.0/(rho_sum*du);
    for (int i = 0; i < nrho; i++) rho[i] *= scal;
}

#endif
