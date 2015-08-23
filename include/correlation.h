#ifndef CORRELATION_H
#define CORRELATION_H

void correlation(int N, const double *u, int ncorr, double *corr);
void meanSqrDrift(int N, const double *u, int ncorr, double *du2);

void correlation_direct(int N, const double *u, int ncorr, double *corr);
void meanSqrDrift_direct(int N, const double *u, int ncorr, double *du2);

#endif
