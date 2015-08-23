#include "cxxheaders.h"
#include "shearsys.h"
#include "geom_oper.h"
#include "ewald.h"
#include "mblas.h"
#include "param.h"

/* ShearSys::addCollisionForce */


namespace {
    struct Cld {
	// imesh -- index of the mesh on which the collision force is applied
	// jmesh -- index of the other mesh in collision
	// x -- contact point
	// nrml -- directional vector (unit length, and points from jmesh to imesh)
	// dist -- separation distance
        int imesh, jmesh;
	double x[3], nrml[3], dist;
	int next;

	Cld() : imesh(-1), jmesh(-1), dist(FLT_MAX), next(-1)
	{
	    x[0] = x[1] = x[2] = 0.0;
	    nrml[0] = nrml[1] = nrml[2] = 0.0;
	}

    };

    // Collision force magnitude
    // Arguments:
    //  x = distance normalized by threshhold
    double CldForceMag(double x) {
        x= fabs(x);
	return (x < 1.0) ? 1.0/(exp(x) - 1.0) - 1.0/(exp(1.0) - 1.0) : 0.0 ;
    }

    // Compare two collision
    bool CldCmp(const Cld &c1, const Cld &c2) {
        return (c1.imesh < c2.imesh) || 
	       (c1.imesh == c2.imesh && c1.jmesh < c2.jmesh) ||
	       (c1.imesh == c2.imesh && c1.jmesh == c2.jmesh && c1.dist < c2.dist);
    }

    void updateCldInfo(NbrList &nlist, double DIST_EPS, MArray<int,1> &firstCld, vector<Cld> &clds);
};


/* Update mesh-mesh collisioninformation
 * Arguments:
 *   nlist -- neighbor list
 *   DIST_EPS -- separation distance threshhold
 *   firstCld -- firstCld[i] points to the first collision record in clds[:] for
 *   		the i-th mesh
 *   clds -- list of collisions
 */
void updateCldInfo(NbrList &nlist, double DIST_EPS, MArray<int,1> &firstCld, vector<Cld> &clds)
{
    for (int foo1 = 0; foo1 < nlist.verts.size(); foo1++) {
	for (int foo2 = nlist.firstNbr[foo1]; foo2 < nlist.firstNbr[foo1+1]; foo2++) {
	    Point *vert = nlist.verts[foo1];
	    Tri *face = nlist.faces[foo2];

	    // exclude self-collision
	    if (vert->mesh == face->mesh) continue;
	    if (nlist.dists[foo2] > DIST_EPS) continue;

	    // Calculate the seperation vector
	    double xtar[3], xtri[3][3];
	    double dist, s0, t0;

	    FOR_I3 xtar[i] = vert->x[i];
	    FOR_I3 FOR_J3 xtri[i][j] = face->vert[i]->x[j];
	    ewald::to_CloseBy(xtri[0], xtar);
	    dist = minDistToTri(xtar, xtri, s0, t0);
	    if (dist > DIST_EPS) continue;

	    double nrml[3];
	    FOR_I3 {
		nrml[i] = xtar[i] - ((1-s0-t0)*xtri[0][i] + s0*xtri[1][i] + t0*xtri[2][i]);
	    }
	    normalizeVec3D(nrml);

	    // Update the collision list
	    int imesh = vert->mesh->Gindx;
	    int jmesh = face->mesh->Gindx;

	    int p = firstCld(imesh);

	    while (p >= 0) {
		if (clds[p].jmesh == jmesh)
		    break;
		else
		    p = clds[p].next;
	    }

	    if (p >= 0) {
		Cld &cld = clds[p];

		if (dist < cld.dist) {
		    FOR_I3 cld.x[i] = vert->x[i];
		    FOR_I3 cld.nrml[i] = nrml[i];
		    cld.dist = dist;
		}
	    } else if (dist < DIST_EPS) {
		Cld newcld;
		newcld.imesh = imesh;
		newcld.jmesh = jmesh;
		newcld.dist = dist;
		FOR_I3 newcld.x[i] = vert->x[i];
		FOR_I3 newcld.nrml[i] = nrml[i];
		newcld.next = firstCld(imesh);

		clds.push_back(newcld);
		firstCld(imesh) = clds.size() - 1;
	    }
	} // foo2
    } // foo1
}


/* Calculate the extra force and torque on particles due to collision
 * Arguments:
 *  iblk -- the mesh block where collision forces are applied
 *  jblk -- the object of collision
 *  DIST_EPS -- distance thresh hold
 * Note:
 *  - The collision force is added to the f[] array of every mesh points */
void ShearSys::addCollisionForce(int iblk, int jblk, double DIST_EPS)
{
    int mpi_rank, mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // Init
    int nmesh = meshList[iblk].size();
    double (*F)[3] = new double[nmesh][3];
    double (*T)[3] = new double[nmesh][3];

    std::fill_n(*F, 3*nmesh, 0.0);
    std::fill_n(*T, 3*nmesh, 0.0);

    // Each processor finds the local collision information
    vector<Cld> myclds;
    // For the i-th mesh, let
    //   p0 = firstCld[i]
    //   p1 = clds[p0].next
    //   p2 = clds[p1].next
    //   ....
    //   and the collision information is stored by clds[p0], clds[p1], clds[p2], ...
    MArray<int,1> firstCld(nmesh);
    firstCld = -1;

    if (ewald::phys::active) {
        NbrList &nlist = nlist_phys[iblk][jblk];
        ::updateCldInfo(nlist, DIST_EPS, firstCld, myclds);
    }

    // Root gathers all collision info
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

    // Root calculates collision force
    if (mpi_rank == 0) {
	// First sort the collision information
	std::sort(recvbuf, recvbuf+displs[mpi_size], CldCmp);

	double dist_min = FLT_MAX;
	int imesh_min = -1;
	int jmesh_min = -1;

	int last_imesh = -1;
	int last_jmesh = -1;
	for (int p = 0; p < displs[mpi_size]; p++) {
	    Cld &cld = recvbuf[p];
	    int imesh = cld.imesh;
	    int jmesh = cld.jmesh;

	    // Each collision is compared only once
	    if (imesh == last_imesh && jmesh == last_jmesh)  {
	        continue;
	    } 
	    else {
	        last_imesh = imesh;
		last_jmesh = jmesh;
            }

	    // The collision force is along nrml[]
	    // The hydrodynamic force is then along -nrml[] for counter balance
	    double f = 5*CldForceMag(cld.dist/DIST_EPS);

	    Mesh *mesh = meshList[iblk][imesh];
	    double xx[3], dF[3], dT[3];
	    FOR_I3 {
	        xx[i] = cld.x[i] - mesh->center[i];
	        dF[i] = f*cld.nrml[i];
	    }
	    cross_product(xx, dF, dT);

	    FOR_I3 F[imesh][i] -= dF[i];
	    FOR_I3 T[imesh][i] -= dT[i];

	    if (cld.dist < dist_min) {
	        dist_min = cld.dist;
		imesh_min = imesh;
		jmesh_min = jmesh;
	    }
	}

	if (dist_min < DIST_EPS) {
	    printf("    %1d -- %1d min sep: %3d %3d %12.2E\n", 
	    		iblk, jblk, imesh_min, jmesh_min, dist_min);
	}
    }

    // Broadcast collision force
    MPI_Bcast(*F, 3*nmesh, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(*T, 3*nmesh, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* Apply collision force
     * For every mesh point, f = ft + fr x (x - xc)
     *  -- ft is a uniform surface force density
     *   -- fr is for torque
     * F = area * ft
     * T = M * fr, where M is the momentum of inertia tensor */
    for (int imesh = 0; imesh < nmesh; imesh++) {
        Mesh *mesh = meshList[iblk][imesh];

	double ft[3], fr[3];
	FOR_I3 ft[i] = F[imesh][i]/mesh->area;
	FOR_I3 fr[i] = mesh->paxis[0][i]*m_ddot(3, mesh->paxis[0], T[imesh])/mesh->pmom[0]
	             + mesh->paxis[1][i]*m_ddot(3, mesh->paxis[1], T[imesh])/mesh->pmom[1]
		     + mesh->paxis[2][i]*m_ddot(3, mesh->paxis[2], T[imesh])/mesh->pmom[2];

        for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    Point &vert = mesh->verts[ivert];

	    // Force
	    FOR_I3 vert.f[i] += ft[i];

	    // Torque
	    double xx[3], df[3];
	    FOR_I3 xx[i] = vert.x[i] - mesh->center[i];
	    cross_product(fr, xx, df);
	    FOR_I3 vert.f[i] += df[i];
	}
    }

    // delete temp arrays
    delete [] F;
    delete [] T;
    delete [] sendbuf;
    delete [] recvbuf;
    delete [] nrecv;
    delete [] displs;
}
