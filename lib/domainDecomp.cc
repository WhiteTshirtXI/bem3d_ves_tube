#include "cxxheaders.h"
#include "domainDecomp.h"

/* domainDecomp_ORB */

namespace {
    // Operator class for comparing between two IntPoint objects
    class IntPointCmp {
        public:
	    int _d;		// the direction for comparison
	    IntPointCmp(); 	// disallow default constructor
	    IntPointCmp(int d) : _d(d) {}
	    ~IntPointCmp() {}

	    bool operator()(IntPoint p0, IntPoint p1) 
	    { 
	        return p0.x[_d] < p1.x[_d]; 
	    }
    };
};


/* Domain decomposition using orthogonal recursive bisection
 * Arguments:
 *  pts -- points
 *  ndims -- number of dimensions
 *  nprocs -- number of processors
 *  myrank -- rank in processors 
 * Note:
 *  1. pts is of integer type to avoid any possible float point error */
void domainDecomp_ORB(vector<IntPoint> &pts, int ndims, int nprocs, int myrank)
{
    if (pts.size() == 0 || nprocs <= 1) return;

    // Find the longest side to cut
    vector<int> bmin(ndims), bmax(ndims);

    for (int d = 0; d < ndims; d++) {
        bmin[d] = INT_MAX;
	bmax[d] = INT_MIN;
    }

    for (int i = 0; i < pts.size(); i++) {
        for (int d = 0; d < ndims; d++) {
	    int xtmp = pts[i].x[d];
	    if (xtmp < bmin[d]) bmin[d] = xtmp;
	    if (xtmp > bmax[d]) bmax[d] = xtmp;
	}
    }

    int dcut = 0;
    for (int d = 1; d < ndims; d++) {
        if ( bmax[d] - bmin[d] > bmax[dcut] - bmin[dcut] ) dcut = d;
    }

    // Cut the domain size in half
    int nprocs1 = nprocs/2;
    int npts1 = (pts.size()*nprocs1)/nprocs;

    IntPointCmp cmp(dcut);
    std::nth_element(pts.begin(), pts.begin()+npts1, pts.end(), cmp);

    // Split the vector
    if (myrank < nprocs1) {
        pts.resize(npts1);
        domainDecomp_ORB(pts, ndims, nprocs1, myrank);
    }
    else {
        pts.erase(pts.begin(), pts.begin()+npts1);
        domainDecomp_ORB(pts, ndims, nprocs-nprocs1, myrank-nprocs1);
    }
}


/* 1D domain decomposition */
void domainDecomp_ORB(int N, const int *x, vector<int> &idx, int nprocs, int myrank)
{
    vector<IntPoint> pts(N);
    for (int i = 0; i < N; i++) {
        pts[i].x[0] = x[i];
	pts[i].idx = i;
    }

    domainDecomp_ORB(pts, 1, nprocs, myrank);

    int nloc = pts.size();
    idx.resize(nloc);
    for (int i = 0; i < nloc; i++) idx[i] = pts[i].idx;
}


/* 2D domain decomposition */
void domainDecomp_ORB(int N, const int (*x)[2], vector<int> &idx, int nprocs, int myrank)
{
    vector<IntPoint> pts(N);
    for (int i = 0; i < N; i++) {
        pts[i].x[0] = x[i][0];
        pts[i].x[1] = x[i][1];
	pts[i].idx = i;
    }

    domainDecomp_ORB(pts, 2, nprocs, myrank);

    int nloc = pts.size();
    idx.resize(nloc);
    for (int i = 0; i < nloc; i++) idx[i] = pts[i].idx;
}


/* 3D domain decomposition */
void domainDecomp_ORB(int N, const int (*x)[3], vector<int> &idx, int nprocs, int myrank)
{
    vector<IntPoint> pts(N);
    for (int i = 0; i < N; i++) {
        pts[i].x[0] = x[i][0];
        pts[i].x[1] = x[i][1];
        pts[i].x[2] = x[i][2];
	pts[i].idx = i;
    }

    domainDecomp_ORB(pts, 3, nprocs, myrank);

    int nloc = pts.size();
    idx.resize(nloc);
    for (int i = 0; i < nloc; i++) idx[i] = pts[i].idx;
}
