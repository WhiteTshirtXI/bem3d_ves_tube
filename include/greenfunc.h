#ifndef GREENFUNC_H
#define GREENFUNC_H

#include "tri.h"

typedef void (*pGreenFunc)(const double *x0, const double *x, const double *n,
		double (*G)[3], const double *f, double *fG, 
		double (*T)[3], const double *g, double *gT);
		
void greenf_unbound(const double *x0, const double *x, const double *n,
		double (*G)[3], const double *f, double *fG, 
		double (*T)[3], const double *g, double *gT);

void greenf_plate(const double *x0, const double *x, const double *n,
		double (*G)[3], const double *f, double *fG, 
		double (*T)[3], const double *g, double *gT);


void greenf_int_on_tri_regular(pGreenFunc, const double *x0, const Tri&,
		double, double*, double(*)[3][3],
		double, double*, double(*)[3][3],
            	const char *qrule="TRI_3");

void greenf_int_on_tri_singular(pGreenFunc, const double *x0, const Tri&,
		double, double*, double(*)[3][3],
		double, double*, double(*)[3][3],
		int nq=4);

#endif
