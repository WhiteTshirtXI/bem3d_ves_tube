/* Source list
 * Note:
 *   -- For rectangular domain, if the side length of the bounding box 
 *      is bigger than L - 2*rc, we make it periodic in that direction
 *   -- For non-orthogonal domain, the bound box is never periodic */

#include "cxxheaders.h"
#include "point.h"
#include "tri.h"
#include "ewald.h"
#include "mathfunc.h"
#include "domainDecomp.h"

/* SrcList::build
 * SrcList::buildBlks
 * SrcList::buildBlks_rect
 * SrcList::buildBlks_general
 * SrcList::getBlkNum 
 * SrcList::getBlkNum_rect 
 * SrcList::getBlkNum_general */

using ewald::L;
using ewald::iL;
using ewald::phys::rc;

using ewald::bctype;
using ewald::RECT;
using ewald::SHEAR;
using ewald::STRAIN;

/* Domain decomposition
 * Arguments:
 *   allFaces -- global faces
 * Note:
 *   -- Local faces are stored in faces[] */
void SrcList::build(const vector<Tri*> &allFaces)
{
    assert(ewald::phys::active);

    using ewald::phys::comm;
    int mpi_size, mpi_rank;
    MPI_Comm_size(comm, &mpi_size);
    MPI_Comm_rank(comm, &mpi_rank);

    // Hashing
    const int IMAX = 2048;
    int N = allFaces.size();
    int (*xint)[3] = new int[N][3];
    for (int i = 0; i < N; i++) {
	FOR_D3 {
	    double xtmp = THRD*(allFaces[i]->vert[0]->x[d]
		    + allFaces[i]->vert[1]->x[d]
		    + allFaces[i]->vert[2]->x[d] );
	    xtmp *= iL[d];
	    xint[i][d] = (int)nearbyint(IMAX*xtmp);
	}
    }
    MPI_Bcast(*xint, 3*N, MPI_INT, 0, comm);

    // Domain decomposition
    vector<int> idx;
    domainDecomp_ORB(N, xint, idx, mpi_size, mpi_rank);
    delete [] xint;

    faces.resize(idx.size());
    for (int i = 0; i < idx.size(); i++) faces[i] = allFaces[idx[i]];

    buildBlks();
}


/* Build Verlet cell list */
void SrcList::buildBlks()
{
    if (bctype == RECT) 
	buildBlks_rect();
    else
	buildBlks_general();
}


/* Build blocks for rectangular domain */
void SrcList::buildBlks_rect()
{
    // Trivial case
    if (numFaces() == 0) {
        FOR_D3 lb[d] = 0.0;
        FOR_D3 ub[d] = 0.0;
	FOR_D3 nblk[d] = 0;

	next.resize(0);
	hoc.resize(0,0,0);
	return;
    }

    // Bounding box
    FOR_D3 lb[d] = FLT_MAX;
    FOR_D3 ub[d] = -FLT_MAX;

    for (int i = 0; i < numFaces(); i++) {
	for (int j = 0; j < 3; j++) {
	    FOR_D3 {
		double xtmp = faces[i]->vert[j]->x[d];
		lb[d] = min(xtmp, lb[d]);
		ub[d] = max(xtmp, ub[d]);
	    }
	}
    } 

    // Extend the box to avoid float point error
    FOR_D3 {
        double eps = 1.E-4*L[d];
        lb[d] -= eps;
        ub[d] += eps;
    }

    // The box size must be at least rc
    FOR_D3 {
        if (ub[d] - lb[d] < rc) {
	    double xx = rc - (ub[d] - lb[d]);
	    lb[d] -= 0.5*xx;
	    ub[d] += 0.5*xx;
	}
    }

    // If box size is greater than (L-2*rc), it becomes periodic in that direction
    FOR_D3 {
	if (ub[d] - lb[d] >= L[d] - 2*rc) {
	    lb[d] = 0.0;
	    ub[d] = L[d];
	    cyclic[d] = true;
	} else {
	    cyclic[d] = false;
        }
    }

    // Create blocks
    FOR_D3 {
	nblk[d] = (int)floor((ub[d] - lb[d])/rc);
	nblk[d] = max(nblk[d], 1);

	blkSz[d] = (ub[d] - lb[d])/nblk[d];
	iblkSz[d] = 1.0/blkSz[d];
    }


    // Hashing
    next.resize(numFaces()); 

    {
	IndexRange r0(-1,nblk[0]+1);
	IndexRange r1(-1,nblk[1]+1);
	IndexRange r2(-1,nblk[2]+1);
	hoc.resize(r0, r1, r2);
    }

    next = -1;
    hoc = -1;

    for (int i = 0; i < numFaces(); i++) {
	Tri *T = faces[i];

	double xc[3];
	FOR_D3 xc[d] = THRD*(T->vert[0]->x[d]
		           + T->vert[1]->x[d]
		           + T->vert[2]->x[d]);

	int i0, i1, i2;
	getBlkNum(xc, i0, i1, i2);
	assert (i0 >= 0 && i0 < nblk[0] &&
		i1 >= 0 && i1 < nblk[1] &&
		i2 >= 0 && i2 < nblk[2]);
	
//	int i0, i1, i2;
//	getBlkNum(xc,i0,i1,i2);
//	assert (i0 < nblk[0]);
//	assert (i1 < nblk[1]);
//	assert (i2 < nblk[2]);
//	assert (i0 >= 0);
//	assert (i1 >= 0);
	assert (i2 >= 0);

	next(i) = hoc(i0,i1,i2);
	hoc(i0,i1,i2) = i;
    }

    // Periodicity
    if (cyclic[0] && nblk[0] > 2) {
	for (int i1 = -1; i1 <= nblk[1]; i1++)
	    for (int i2 = -1; i2 <= nblk[2]; i2++) {
		hoc(-1,i1,i2) = hoc(-1+nblk[0],i1,i2);
		hoc(nblk[0],i1,i2) = hoc(0,i1,i2);
	    }
    }

    if (cyclic[1] && nblk[1] > 2) {
	for (int i0 = -1; i0 <= nblk[0]; i0++)
	    for (int i2 = -1; i2 <= nblk[2]; i2++) {
		hoc(i0,-1,i2) = hoc(i0,-1+nblk[1],i2);
		hoc(i0,nblk[1],i2) = hoc(i0,0,i2);
	    }
    }


    if (cyclic[2] && nblk[2] > 2) {
	for (int i0 = -1; i0 <= nblk[0]; i0++)
	    for (int i1 = -1; i1 <= nblk[1]; i1++) {
		hoc(i0,i1,-1) = hoc(i0,i1,-1+nblk[2]);
		hoc(i0,i1,nblk[2]) = hoc(i0,i1,0);
	    }
    }
}


void SrcList::buildBlks_general()
{
    // Trivial case
    if (numFaces() == 0) {
        FOR_D3 lb[d] = 0.0;
        FOR_D3 ub[d] = 0.0;
	FOR_D3 nblk[d] = 0;

	next.resize(0);
	hoc.resize(0,0,0);
	return;
    }

    // Bounding box
    FOR_D3 lb[d] = FLT_MAX;
    FOR_D3 ub[d] = -FLT_MAX;

    for (int i = 0; i < numFaces(); i++) {
	for (int j = 0; j < 3; j++) {
	    FOR_D3 {
		double xtmp = faces[i]->vert[j]->x[d];
		lb[d] = min(xtmp, lb[d]);
		ub[d] = max(xtmp, ub[d]);
	    }
	}
    } 

    // Extend the box size to avoid float point error
    FOR_D3 {
        double eps = 1.E-4*L[d];
	lb[d] -= eps;
	ub[d] += eps;
    }

    // Box size must be at least rc
    FOR_D3 {
        if (ub[d] - lb[d] < rc) {
	    double xx = rc - (ub[d] - lb[d]);
	    lb[d] -= 0.5*xx;
	    ub[d] += 0.5*xx;
	}
    }

    // Create blocks
    FOR_D3 {
	nblk[d] = (int) floor((ub[d] - lb[d])/rc);
	nblk[d] = max(nblk[d], 1);

	blkSz[d] = (ub[d] - lb[d])/nblk[d];
	iblkSz[d] = 1.0/blkSz[d];
    }


    // Hashing
    next.resize(numFaces());

    {
	IndexRange r0(0,nblk[0]);
	IndexRange r1(0,nblk[1]);
	IndexRange r2(0,nblk[2]);
	hoc.resize(r0, r1, r2);
    }

    next = -1;
    hoc = -1;

    for (int i = 0; i < numFaces(); i++) {
	Tri *T = faces[i];

	double xc[3];
	FOR_D3 xc[d] = THRD*(T->vert[0]->x[d]
	                   + T->vert[1]->x[d]
	                   + T->vert[2]->x[d]);

	int i0, i1, i2;
	getBlkNum(xc, i0, i1, i2);
	assert (i0 >= 0 && i0 < nblk[0] &&
		i1 >= 0 && i1 < nblk[1] &&
		i2 >= 0 && i2 < nblk[2]);

	next(i) = hoc(i0,i1,i2);
	hoc(i0,i1,i2) = i;
    }
}


/* Compute the cell block number 
 * Arguments:
 *   x -- coordinate
 *   i0,i1,i2 -- block number */
void SrcList::getBlkNum(const double *x, int &i0, int &i1, int &i2)
{
    if (ewald::bctype == ewald::RECT)
	getBlkNum_rect(x, i0, i1, i2);
    else
	getBlkNum_general(x, i0, i1, i2);
}


void SrcList::getBlkNum_rect(const double *x, int &i0, int &i1, int &i2)
{
    int ix[3];

    FOR_D3 {
	if (cyclic[d]) {
	    ix[d] = (int) floor(x[d]*iblkSz[d]);
	    ix[d] = modulo(ix[d], nblk[d]);
	}
	else {
	    double xx = x[d] - 0.5*(lb[d] + ub[d]);
	    xx = nearbyint(xx*iL[d])*L[d];
	    ix[d] = (int) floor((x[d] - xx - lb[d])*iblkSz[d]);
	}
    }

    i0 = ix[0];
    i1 = ix[1];
    i2 = ix[2];
}


void SrcList::getBlkNum_general(const double *x, int &i0, int &i1, int &i2)
{
    i0 = (int) floor( (x[0] - lb[0])*iblkSz[0] );
    i1 = (int) floor( (x[1] - lb[1])*iblkSz[1] );
    i2 = (int) floor( (x[2] - lb[2])*iblkSz[2] );
}
