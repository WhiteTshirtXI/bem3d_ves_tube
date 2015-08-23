// PETSc wrapper

#ifndef MYPETSC_H
#define MYPETSC_H

#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"
#include "petscpc.h"

void myAllReduce(int nvar, double *u, int nloc, double *uloc, int *l2g, MPI_Comm comm);

void myVecSetValues(Vec, const double *, InsertMode);
void myVecScatterCreateToAll(Vec, VecScatter&);
void myVecScatter(VecScatter, Vec, Vec, InsertMode, ScatterMode);
void myVecScatter(VecScatter, Vec, double*, InsertMode, ScatterMode);
void myVecScatter(VecScatter, double*, Vec, InsertMode, ScatterMode);

void myMatMult(Mat, double *, double *);
void myMatMultAdd(Mat, double *, double *);
void myMatMultAdd(double, Mat, double *, double *);

#endif
