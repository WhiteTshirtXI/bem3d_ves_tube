#ifndef MEMBRANE_H
#define MEMBRANE_H

#include "cxxheaders.h"
#include "classdec.h"

void tri_strain_energy(Tri &, Tri &, double, double, double *, double (*)[3]);
void edge_bend_energy(Edge &, double, double, double *, double(*)[3], double (*)[3]);

#endif
