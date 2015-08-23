#include "cxxheaders.h"
#include "stokesys.h"
#include "ewald.h"
#include "mblas.h"
#include "param.h"
#include "geom_oper.h"
#include "debugfunc.h"

/* StokeSys::cellNoContact
 * updateCldInfo */


namespace {
    // collision information
    // Use a link-list to store collision info for the same mesh point
    struct Cld {
        int ivert;
        double dist, nrml[3];

	Cld() : ivert(-1), dist(FLT_MAX) 
	{
	    nrml[0] = nrml[1] = nrml[2] = 0.0;
	}

    };

    // Compare two collision
    bool CldCmp(const Cld &c1, const Cld &c2) {
        return ( (c1.ivert < c2.ivert) || 
	         (c1.ivert == c2.ivert && c1.dist < c2.dist) );
    }

    void updateCldInfo(NbrList &nlist, double DIST_EPS, 
    		MArray<int,1> &pcld, vector<Cld> &clds);
};


/* Update collision information
 * Arguments:
 *   nlist -- the neighborlist
 *   DIST_EPS -- threshhold distance 
 *   pcld -- clds[pcld[i]] stores the collision info of the i-th point
 *   clds -- collisions */
void updateCldInfo(NbrList &nlist, double DIST_EPS, MArray<int,1> &pcld, vector<Cld> &clds)
{
    for (int foo1 = 0; foo1 < nlist.verts.size(); foo1++) {
        Point &vert = *nlist.verts[foo1];
	int ivert = vert.Gindx;

	double min_dist = FLT_MAX;
	double min_nrml[3];

        for (int foo2 = nlist.firstNbr[foo1]; foo2 < nlist.firstNbr[foo1+1]; foo2++) {
	    // Find the accurate distance to the face
	    Tri &face = *nlist.faces[foo2];

	    // exclude mesh self-contact
	    if (vert.mesh != NULL && vert.mesh == face.mesh) continue;	
	    if (nlist.dists[foo2] > DIST_EPS || nlist.dists[foo2] > min_dist) continue;
	    
	    double xi[3];
	    m_dcopy(3, vert.x, xi);
	    ewald::to_CloseBy(face.xc, xi);

	    double xtri[3][3];
	    for (int l = 0 ; l < 3; l++) m_dcopy(3, face.vert[l]->x, xtri[l]);

	    double rr, s0, t0;
	    rr = minDistToTri(xi, xtri, s0, t0);

	    // Update minimum separation and outward normal directions
	    if (rr < min_dist) {
	        min_dist = rr;

		FOR_D3 {
		    double xtmp = (1.0-s0-t0)*xtri[0][d] + s0*xtri[1][d] + t0*xtri[2][d];
		    min_nrml[d] = xi[d] - xtmp;
		}
		normalizeVec3D(min_nrml);
	    }
	} // foo2


        if (min_dist > DIST_EPS) continue;


	// Update collision list
	int p = pcld(ivert);
	if (p >= 0) {
	    Cld &cld = clds[p];

	    if (min_dist < cld.dist) {
		cld.dist = min_dist;
		FOR_I3 cld.nrml[i] = min_nrml[i];
	    }
	} 
	else {
	    Cld newcld;
	    newcld.ivert = ivert;
	    newcld.dist = min_dist;
	    FOR_D3 newcld.nrml[d] = min_nrml[d];

	    clds.push_back(newcld);
	    pcld(ivert) = clds.size() - 1;
        }
    } // ivert
}


/* Prevent cells from contacting each other
 * Arguments:
 *   DIST_EPS -- distance thresh hold */
void StokeSys::cellNoContact(double DIST_EPS)
{
    int mpi_rank, mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // Get local collision information
    int nvert = vertList[CELL].size();
    vector<Cld> myclds;
    MArray<int,1> pcld(nvert);
    pcld = -1;

    if (ewald::phys::active) {
	::updateCldInfo(nlist_phys[CELL][CELL], DIST_EPS, pcld, myclds);
	::updateCldInfo(nlist_phys[CELL][RIGID], DIST_EPS, pcld, myclds);
	::updateCldInfo(nlist_phys[CELL][WALL], DIST_EPS, pcld, myclds);
    }

    // root collect all collision information
    int nsend = myclds.size();
    int *nrecv = new int[mpi_size]; 
    int *displs = new int[mpi_size+1];
    MPI_Gather(&nsend, 1, MPI_INT, nrecv, 1, MPI_INT, 0, MPI_COMM_WORLD);
    displs[0] = 0;
    for (int i = 1; i <= mpi_size; i++) {
	displs[i] = displs[i-1] + nrecv[i-1];
    }

    Cld *sendbuf = new Cld[nsend];
    std::copy(myclds.begin(), myclds.end(), sendbuf);

    Cld *recvbuf = NULL;
    if (mpi_rank == 0) recvbuf = new Cld[displs[mpi_size]];

    MPI_Datatype MPI_CLD;
    MPI_Type_contiguous(sizeof(Cld), MPI_CHAR, &MPI_CLD);
    MPI_Type_commit(&MPI_CLD);

    MPI_Gatherv(sendbuf, nsend, MPI_CLD,
    		recvbuf, nrecv, displs, MPI_CLD, 
		0, MPI_COMM_WORLD);
    MPI_Type_free(&MPI_CLD);

    double (*dx)[3] = new double[nvert][3];
    m_dclear(3*nvert, *dx);

    if (mpi_rank == 0) {
	// First sort the collision information
	std::sort(recvbuf, recvbuf+displs[mpi_size], CldCmp);

	double dist_min = FLT_MAX;
	int cnt = 0;

	int last_ivert = -1;
	for (int foo = 0; foo < displs[mpi_size]; foo++) {
	    Cld &cld = recvbuf[foo];
	    int ivert = cld.ivert;

	    // Each point is moved only once using the smallest separation
	    if (ivert == last_ivert)  {
	        continue;
	    } 
	    else {
	        cnt++;
	        last_ivert = ivert;
            }

	    double ds = DIST_EPS - cld.dist;
	    FOR_J3 dx[ivert][j] = ds*cld.nrml[j];

	    dist_min = std::min(dist_min, cld.dist);
	}


	if (cnt > 0) {
	    printf("    cell-cell min sep = %.2E", dist_min);
	    printf("    %d mesh points moved\n", cnt);
	}
    }

    MPI_Bcast(*dx, 3*nvert, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i = 0; i < nvert; i++) {
        Point *vert = vertList[CELL][i];
	m_dadd(3, dx[i], vert->x);
    }

    delete [] nrecv;
    delete [] displs;
    delete [] sendbuf;
    delete [] recvbuf;
    delete [] dx;
}


/* Check whether a point lies within a cell 
 * Argument:
 *   x0 -- the coordinate of the point
 *   whichCell -- index of the cell that the point penetrates 
 * Algorithm:
 *   -- If a point is outside a cell, the total spherical angle is 0 */
bool StokeSys::pointInsideSomeCell(const double *x0, int *whichCell)
{
    bool is_interior = false;
    if (whichCell) *whichCell = -1;

    for (int icell = 0; icell < numCells(); icell++) {
	Vesicle &cell = cells[icell];

	double xtmp[3];
	FOR_I3 xtmp[i] = x0[i];
	ewald::to_CloseBy(cell.center, xtmp);

	double xmin[3], xmax[3];
	cell.getCoordRange(xmin, xmax);

	if (   xtmp[0] < xmin[0] || xtmp[0] > xmax[0] 
	    || xtmp[1] < xmin[1] || xtmp[1] > xmax[1]
	    || xtmp[2] < xmin[2] || xtmp[2] > xmax[2] ) continue;

	double sangle = 0.0;
	for (int iface = 0; iface < cell.numFaces(); iface++) {
	    Tri &face = cell.faces[iface];
	    double xtri[3][3];
	    FOR_I3 FOR_J3 xtri[i][j] = face.vert[i]->x[j];
	    sangle += tri_solidAngle(xtmp, xtri);
	}
	if (fabs(sangle) > 1.E-2) {
	    is_interior = true;
	    if (whichCell) *whichCell = icell;
	    break;
	}
    } // icell

    return is_interior;
}
