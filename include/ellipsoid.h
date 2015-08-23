#ifndef ELLIPSOID_H
#define ELLIPSOID_H

namespace ellipsoid 
{
    void calcLatLonAngles(double A, double B, double C, const double *x,
		double &th, double &phi);

    void calcLocalGeom(double A, double B, double C, double th, double phi,
    		double *x, double *nrml, double &H, double &K);
};

#endif
