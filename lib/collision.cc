#include "cxxheaders.h"
#include "collision.h"
#include "mesh.h"
#include "mathfunc.h"

/* collision::updateCldInfo 
 * collision::forceSeparation
 */

/* Update collision information
 * Arguments:
 *   nlist -- the neighborlist
 *   DIST_EPS -- threshhold distance 
 *   pcld -- clds[pcld[i]] stores the collision info of the i-th point
 *   clds -- collisions 
 * Note:
 *   -- For each mesh, we define a shell surface that is a small distance away from the 
 *      mesh surface, the points within the shell are moved onto the surface.
 */
void collision::updateCldInfo(NbrList &nlist, double DIST_EPS, MArray<int,1> &pcld, vector<Cld> &clds)
{
    for (int foo1 = 0; foo1 < nlist.verts.size(); foo1++) {
        Point &vert = *nlist.verts[foo1];
	int ivert = vert.Gindx;

	double min_spr = FLT_MAX;
	double dx[3];

        for (int foo2 = nlist.firstNbr[foo1]; foo2 < nlist.firstNbr[foo1+1]; foo2++) {
	    // Find the accurate distance to the face
	    Tri &face = *nlist.faces[foo2];

	    // Exclude mesh self-contact
	    if (vert.mesh != NULL && vert.mesh == face.mesh) continue;	

	    // Find if the point is within the shell surface
	    if (nlist.dists[foo2] > DIST_EPS) continue;

	    double xtri[3][3], xi[3], s0, t0, spr;

	    for (int l = 0; l < 3; l++) {
	        int ivert = face.ivert[l];
		m_dcopy(3, face.vert[l]->x, xtri[l]);
            }

	    m_dcopy(3, vert.x, xi);
	    ewald::to_CloseBy(xtri[0], xi);

	    spr = minDistToTri(xi, xtri, s0, t0);
	    if (spr > DIST_EPS || spr > min_spr) continue;

	    // Build an triangle on the virtual shell surface, and find the projection
	    // of xi on the triangle
	    min_spr = spr;

	    for (int l = 0; l < 3; l++) {
	        int ivert = face.ivert[l];
		m_daxpy(3, spr + 0.60*(DIST_EPS - spr), &face.mesh->vertNrml(ivert,0), xtri[l]);
	    }
	    
	    minDistToTri(xi, xtri, s0, t0);

	    FOR_D3 {
		double xtmp = (1.0-s0-t0)*xtri[0][d] + s0*xtri[1][d] + t0*xtri[2][d];
		dx[d] = xtmp - xi[d];
	    }
	} // foo2


        if (min_spr > DIST_EPS) continue;

	// Update collision list
	int p = pcld(ivert);
	if (p >= 0) {
	    Cld &cld = clds[p];

	    if (min_spr < cld.spr) {
		cld.spr = min_spr;
		FOR_I3 cld.dx[i] = dx[i];
	    }
	} 
	else {
	    Cld newcld;
	    newcld.ivert = ivert;
	    newcld.spr = min_spr;
	    FOR_D3 newcld.dx[d] = dx[d];

	    clds.push_back(newcld);
	    pcld(ivert) = clds.size() - 1;
        }
    } // ivert
}


/* Prevent meshes from contacting each other
 * Arguments:
 *   vlist -- all the vertics that can be moved to avoid collision
 *   nlists -- the neighborlists of vertex-and-face pairs
 *   DIST_EPS -- distance threshhold */
void collision::forceSeparation(vector<Point*> &vlist, vector<NbrList*> &nlists, double DIST_EPS)
{
    int mpi_rank, mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // Get local collision information
    int nvert = vlist.size();
    vector<Cld> myclds;
    MArray<int,1> pcld(nvert);
    pcld = -1;

    for (int i = 0; i < nlists.size(); i++) {
	updateCldInfo(*nlists[i], DIST_EPS, pcld, myclds);
    }

    // root collects all collision information
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

	double spr_min = FLT_MAX;
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

	    FOR_J3 dx[ivert][j] = cld.dx[j];

	    spr_min = std::min(spr_min, cld.spr);
	}


	if (cnt > 0) {
	    printf("    cell-cell min sep = %.2E", spr_min);
	    printf("    %d mesh points moved\n", cnt);
	}
    }

    MPI_Bcast(*dx, 3*nvert, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    for (int i = 0; i < nvert; i++) {
	m_dadd(3, dx[i], vlist[i]->x);
    }

    delete [] nrecv;
    delete [] displs;
    delete [] sendbuf;
    delete [] recvbuf;
    delete [] dx;
}
