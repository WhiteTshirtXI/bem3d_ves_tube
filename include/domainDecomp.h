#ifndef DOMAINDECOMP_H
#define DOMAINDECOMP_H

// Points with integer coordinate and associated with an index
// Mainly for sorting and domain decomposition
struct IntPoint
{
    int x[3], idx;
};

void domainDecomp_ORB(vector<IntPoint> &pts, int ndims, int nprocs, int myrank);
void domainDecomp_ORB(int N, const int *x,      vector<int> &idx, int nprocs, int myrank);
void domainDecomp_ORB(int N, const int (*x)[2], vector<int> &idx, int nprocs, int myrank);
void domainDecomp_ORB(int N, const int (*x)[3], vector<int> &idx, int nprocs, int myrank);

#endif
