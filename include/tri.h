#ifndef TRI_H
#define TRI_H

#include "classdec.h"
#include "point.h"

struct Tri 
{
    int ivert[3];
    Point *vert[3];
    Mesh *mesh;

    int Gindx;	// global index

    // xc -- centroid
    // normal -- surface normal
    // detJ = determinant of Jacobian
    // area = 0.5*detJ, area
    // diam = min(the diameter of the circumcircle, the longest edge length)
    double xc[3], normal[3], detJ, area, diam;

    // a -- metric tensor
    // detA -- det(a)
    // arcp -- reciprocal metric tensor
    double a[2][2], detA;
    double arcp[2][2];

    // For compute the approximate distance to a point
    double _normal[4][3], _offset[4];

    Tri();
    Tri(int, int, int);
    Tri(const Tri &);
    ~Tri();

    void updateGeometry();

    void baryCoord(const double *xtar, double &s, double &t) const;
    double minDistToPoint_Approx(const double*) const;
    void gradNormal(double (*)[3][3]) const;
    void gradArea(double (*)[3]) const;
    void secondDerivArea(double (*)[3][3][3]) const;

    double minAngle();
};

#endif
