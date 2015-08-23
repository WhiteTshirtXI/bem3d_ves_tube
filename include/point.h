#ifndef POINT_H
#define POINT_H

#include "classdec.h"

struct Point 
{
    double x[3];	// coordinate

    Point *pp;		// the same point on the periodic domain boundary
    Mesh *mesh;		// the mesh
    int indx; 		// index on the mesh
    int Gindx;		// global index

    double f[3], g[3];	// single- and double-layer density
    double C[3][3];	// double-layer jump term

    Point();
    Point(double, double, double); 
    Point(const Point &);
    Point& operator=(const Point &);
    ~Point();
};

#endif
