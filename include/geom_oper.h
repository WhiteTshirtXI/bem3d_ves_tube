// Basic geometric calculations
#ifndef GEOM_OPER_H
#define GEOM_OPER_H

#include "classdec.h"

double normVec3D(const double *);

double normalizeVec3D(double *);

double distCart3D(const double *, const double *);

void cross_product(const double *, const double *, double *);

double triple_product(const double *, const double *, const double *);

void tri_normal(const double *, const double *, const double *, double *, double &);
void tri_normal(const double (*)[3], double *, double &);

double tri_shapeFactor(const double *, const double *, const double *);
double tri_shapeFactor(const double (*)[3]);

double minDistToTri(const double *, const double (*)[3], double &, double &);

bool point_in_seg(double, double, double, double);
bool seg_overlap(double, double, double, double, double);

double tri_solidAngle(const double *, const double (*)[3]);

void solidAngle(const double *, int , MArray<double,2> &, double &, double *);

double bendAngle(const double *, const double *, const double *, const double *);

void calcRotateMatrix(const double (*)[3], const double (*)[3], double (*)[3]);
void calcRotateMatrix(const double *, const double *, double (*)[3]);

void rotate_vector(const double *, double, double *);

void calcPolygonAreaCenter(int, const double (*)[3], double &, double *);

void slice_triangle(const double (*)[3], 
	const double *, const vector<double> &, 
	int &, int *, double (*)[5][3]);

bool tri_box_collide(const double (*xtri)[3], double *xmin, double *xmax);

#endif
