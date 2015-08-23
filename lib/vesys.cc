#include "cxxheaders.h"
#include "vesys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "quadrature.h"
#include "param.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "miscUtils.h"
#include "debugfunc.h"
#include "collision.h"

/* VeSys::VeSys
 * VeSys::~VeSys 
 * VeSys::initialHook 
 * VeSys::domainDecomp
 * VeSys::updateSourceGeometry
 * VeSys::updateNeighborList
 * VeSys::calcBieMat
 * VeSys::rebox
 * VeSys::syncCoord
 * VeSys::calcStress
 * VeSys::writeAll
 * VeSys::writeVesicles
 * VeSys::writeVesicleCenter
 * VeSys::writeStress
 * VeSys::writeRestart
 * VeSys::readRestart
 * VeSys::volFrac
 * VeSys::calcMaxVelPerturb
 * VeSys::solveVel
 * VeSys::timeInt */

VeSys::VeSys()
{ 
    matSL = NULL;
    matDL = NULL;
}


VeSys::~VeSys()
{ }


void VeSys::initialHook()
{
    // Build mesh list and index meshes
    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
        meshList.push_back(&vesicle);

	vesicle.Gindx = ives;
	vesicle.setInternalPointers();
	vesicle.noPeriodicBoundaries();
	vesicle.buildControlPoints();

	int nvert = vesicle.numVerts();
	vesicle.vertArea.resize(nvert);
	vesicle.vertNrml.resize(nvert,3);
	vesicle.vertH.resize(nvert);
	vesicle.sigma.resize(nvert);
	vesicle.v.resize(nvert,3);
    }

    // Build vertex list, vertex indices, and face list
    int vert_cnt = 0;
    int face_cnt = 0;

    for (int imesh = 0; imesh < meshList.size(); imesh++) {
	Mesh *mesh = meshList[imesh];

	// Vertices
	for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    Point &vert = mesh->verts[ivert];
	    if (vert.pp == NULL) {
		vertList.push_back(&vert);
		vert.Gindx = vert_cnt++;
	    }
	}

	// Periodic boundaries
	for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    Point &vert = mesh->verts[ivert];
	    if (vert.pp != NULL) vert.Gindx = vert.pp->Gindx;
	}

	// Faces
	for (int iface = 0; iface < mesh->numFaces(); iface++) {
	    Tri &face = mesh->faces[iface];
	    faceList.push_back(&face);
	    face.Gindx = face_cnt++;
	}
    } // imesh
}


void VeSys::domainDecomp()
{
    // Build source list
    if (ewald::phys::active) {
        slist_phys.build(faceList);
    }

    if (ewald::four::active) {
	ewald::four::setSources(faceList, slist_four);
    }

    // Build target list
    if (ewald::four::active) {
        ewald::four::setTargets(vertList, tlist_four);
    }

    // Set active and private flags of meshes
    {
	// Init
        for (int imesh = 0; imesh < meshList.size(); imesh++) {
	    Mesh *mesh = meshList[imesh];
	    mesh->isActive = 0;
	    mesh->isPrivate = 0;
	}

	// Active flag
	for (int i = 0; i < slist_phys.numFaces(); i++) {
	    Mesh *mesh = slist_phys.faces[i]->mesh;
	    mesh->isActive = 1;
	}

	for (int i = 0; i < slist_four.size(); i++) {
	    Mesh *mesh = slist_four[i]->mesh;
	    mesh->isActive = 1;
	}


	// Private flag
	int mpi_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	int nmesh = meshList.size();
	MArray<int,1> ownerRank(nmesh), ownerRankMax(nmesh);
	for (int imesh = 0; imesh < nmesh; imesh++) {
	    Mesh *mesh = meshList[imesh];
            ownerRank(imesh) = mesh->isActive ? mpi_rank : -1;
	}

        MPI_Allreduce(ownerRank.data(), ownerRankMax.data(), 
			nmesh, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        for (int imesh = 0; imesh < nmesh; imesh++) {
	    Mesh *mesh = meshList[imesh];
	    mesh->isPrivate = (ownerRankMax(imesh) == mpi_rank);
	}
    }
}


/*   Update the geometry of meshes that are active */
void VeSys::updateSourceGeometry()
{
    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	if (vesicle.isActive) vesicle.updateGeometry();
    }
}


// Update neighbor list
void VeSys::updateNeighborList()
{
    if (ewald::phys::active) {
	nlist_phys.build(vertList, slist_phys);
    }
}


// Compute the single- and double-layer matrices
void VeSys::calcBieMat()
{
    if (matSL) MatDestroy(&matSL);
    if (matDL) MatDestroy(&matDL);

    int nrow = 3*nlist_phys.verts.size(); 
    int ncol = 3*vertList.size();

    MatCreate(PETSC_COMM_SELF, &matSL);
    MatSetType(matSL, MATSEQAIJ);
    MatSetSizes(matSL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);

    MatCreate(PETSC_COMM_SELF, &matDL);
    MatSetType(matDL, MATSEQAIJ);
    MatSetSizes(matDL, PETSC_DECIDE, PETSC_DECIDE, nrow, ncol);

    if (ewald::phys::active) {
        ewald::phys::calcSurfIntMat(nlist_phys, matSL, matDL);
    }
}


void VeSys::rebox()
{
    using ewald::L;
    using ewald::iL;

    for (int i = 0; i < meshList.size(); i++) {
        Mesh *mesh = meshList[i];

        // Determine translating distance
        double xmin[3], xmax[3], xc[3];

        mesh->getCoordRange(xmin, xmax);
	FOR_D3 xc[d] = 0.5*(xmin[d] + xmax[d]);

	bool move_mesh = false;
	FOR_D3 {
            double eps = 0.1*L[d];
	    if (xc[d] < -eps || xc[d] > L[d] + eps) move_mesh = true;
	}

	if (move_mesh) {
	    double origin[3], xcnew[3], xx[3];

	    FOR_D3 origin[d] = 0.5*ewald::L[d];
	    FOR_D3 xcnew[d] = xc[d];
	    ewald::to_CloseBy(origin, xcnew);
	    FOR_D3 xx[d] = xcnew[d] - xc[d];

	    for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
		Point &vert = mesh->verts[ivert];
		FOR_D3 vert.x[d] += xx[d];
	    }
	}
    } // i
}

void VeSys::syncCoord()
{
    // Assemble global coordinate arrays
    int npoint = 0;
    for (int i = 0; i < meshList.size(); i++)
        npoint += meshList[i]->numVerts();

    double (*x)[3] = new double[npoint][3];

    for (int i  = 0, p = 0; i < meshList.size(); i++) {
        Mesh *mesh = meshList[i];
	for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    m_dcopy(3, mesh->verts[ivert].x, x[p]);
	    p++;
	}
    }

    // Broadcast from root node
    MPI_Bcast(*x, 3*npoint, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Update vertex coordinates
    for (int i = 0, p = 0; i < meshList.size(); i++) {
        Mesh *mesh = meshList[i];

	for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    m_dcopy(3, x[p], mesh->verts[ivert].x);
	    p++;
	}
    }

    // Dealloc temp arrays
    delete [] x;
}


/* Compute the particle stress
 * Arguments:
 *   tau -- added stress, normalized by (vesicle vol) * (shear rate)  */
void VeSys::calcStress(double (*tau)[3])
{
    // Init
    m_dclear(9, *tau);

    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	int nvert = vesicle.numVerts();

	MArray<double,2> fb(nvert,3), fs(nvert,3);
	vesicle.bendForce(fb);
	vesicle.tensionForce(vesicle.sigma, fs);

	for (int ifa = 0; ifa < vesicle.faces.size(); ifa++) {
	    Tri &face = vesicle.faces[ifa];

	    double xtri[3][3], ftri[3][3], vtri[3][3]; 
	    for (int l = 0; l < 3; l++) {
		int ivert = face.ivert[l];

		FOR_D3 {
		    xtri[l][d] = vesicle.verts[ivert].x[d];
		    ftri[l][d] = fb(ivert,d) + fs(ivert,d);
		    vtri[l][d] = vesicle.v(ivert,d);
	        }
	    }

	    Quad2D &Q = quadrature::select_rule_2d("TRI_3");
	    for (int iq = 0; iq < Q.n(); iq++) {
	        double s = Q.x(iq), t = Q.y(iq);
			double r = 1.0 - s - t;

	        double xq[3], fq[3], vq[3], dA;
		FOR_I3 {
		    xq[i] = r*xtri[0][i] + s*xtri[1][i] + t*xtri[2][i] - vesicle.center[i];
		    fq[i] = r*ftri[0][i] + s*ftri[1][i] + t*ftri[2][i];
		    vq[i] = r*vtri[0][i] + s*vtri[1][i] + t*vtri[2][i];
	        }

		dA = Q.w(iq)*face.detJ;

		FOR_I3
		FOR_J3 {
		    tau[i][j] += 0.5*(fq[i]*xq[j] + fq[j]*xq[i])*dA;
		    tau[i][j] += (viscRat - 1)*(vq[i]*face.normal[j] 
			                      + vq[j]*face.normal[i])*dA;
		}
	    } // iq
	} // ifa
    } // ives

    // normalize by (total vesicle volume)*(shear rate)
    double cvol = 0.0;
    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
        cvol += vesicle.vol;
    }
    m_dscal(9, 1.0/(cvol*shRate), *tau);
} 


void VeSys::writeAll()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    const char FN_FMT[] = "%s%6.6d%s";
    char token[256], fn[256];
    int nout;
    const int LT_CHUNK = 10000;

    // Vesicles
    strcpy(token, "VESICLE_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/vesicle", lt, ".dat");
	if (mpi_rank == 0) writeVesicles(fn);
    }

    // Restart
    strcpy(token, "RESTART_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/restart", lt, ".dat");
	if (mpi_rank == 0) writeRestart(fn);
    }

    // Center of vesicles
    strcpy(token, "CENTER_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
        sprintf(fn, FN_FMT, "D/center", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
        if (mpi_rank == 0) writeVesicleCenter(fn);
    }

    // Stress
    strcpy(token, "STRESS_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	double tau[3][3];
	calcStress(tau);
	sprintf(fn, FN_FMT, "D/stress", Nt0, ".dat");
	if (mpi_rank == 0) writeStress(fn, tau);
    }
}


/* Write vesicles
 * Arguments:
 *   fn -- file name */
void VeSys::writeVesicles(const char *fn)
{
    if (numVesicles() <= 0) return;

    FILE *file = fopen(fn, "w");
    fprintf(file, "variables = x, y, z, u, v, w, sigma\n");

    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	int nvert = vesicle.numVerts();
	int nface = vesicle.numFaces();

	MArray<double,2> x(nvert,3), v(nvert,3);
	MArray<double,1> sigma(nvert);
	MArray<int,2> f2v(nface,3);

	vesicle.getCoords(x);
	vesicle.getConnectivities(f2v);

	if (vesicle.v.size(0) == nvert && vesicle.v.size(1) == 3) {
	    v = vesicle.v;
	} else {
	    v = 0.0;
	}

	if (vesicle.sigma.size(0) == nvert) {
	    sigma = vesicle.sigma;
	} else {
	    sigma = 0.0;
	}

	fprintf(file, "zone N=%d E=%d F=FEPOINT ET=TRIANGLE\n", nvert, nface);
        
	// coordinates
	for (int ivert = 0; ivert < nvert; ivert++) {
	    fprintf(file, " %10.3e %10.3e %10.3e", x(ivert,0), x(ivert,1), x(ivert,2));
	    fprintf(file, " %10.3e %10.3e %10.3e", v(ivert,0), v(ivert,1), v(ivert,2));
	    fprintf(file, " %10.3e\n", sigma(ivert));
        }

	// connectivity
	for (int iface = 0; iface < nface; iface++) {
	    fprintf(file, " %d %d %d\n", f2v(iface,0)+1, f2v(iface,1)+1, f2v(iface,2)+1);
        }
    }

    fclose(file);
}



void VeSys::writeRestart(const char *fn)
{
    hid_t fid;
    hsize_t dims[2];
    char token[256];

    fid = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dims[0] = 3;
    H5LTmake_dataset_double(fid, "EWALD_L", 1, dims, ewald::L);

    dims[0] = 1;
    H5LTmake_dataset_double(fid, "TIME", 1, dims, &time);

    dims[0] = 1;
    H5LTmake_dataset_int(fid, "LT", 1, dims, &lt);

    dims[0] = 1;
    H5LTmake_dataset_double(fid, "STRAIN", 1, dims, &ewald::strain);

    dims[0] = 1;
    H5LTmake_dataset_double(fid, "SHRATE", 1, dims, &shRate);

    int nvesicle = numVesicles();
    if (nvesicle > 0) {
	dims[0] = 1;
	H5LTmake_dataset_int(fid, "NVESICLE", 1, dims, &nvesicle);
    }

    // Vesicles 
    for (int ives = 0; ives < numVesicles(); ives++) {
        Vesicle &vesicle = vesicles[ives];
	int nvert = vesicle.numVerts();
	int nface = vesicle.numFaces();

	MArray<double,2> x(nvert,3);
	MArray<int,2> f2v(nface,3);
	vesicle.getCoords(x);
	vesicle.getConnectivities(f2v);

        sprintf(token, "VESICLE%d", ives);
	hid_t gid = H5Gcreate(fid, token, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	dims[0] = nvert;
	dims[1] = 3;
        H5LTmake_dataset_double(gid, "X", 2, dims, x.data());

	dims[0] = nface;
	dims[1] = 3;
	H5LTmake_dataset_int(gid, "F2V", 2, dims, f2v.data());

	H5Gclose(gid);
    }

    H5Fclose(fid);
}


/* Read restart file in HDF5 format
 * Arguments:
 *  fn -- file name */
void VeSys::readRestart(const char *fn)
{
    const MPI_Comm mpi_comm = MPI_COMM_WORLD;
    int mpi_size = 1;
    int mpi_rank = 0;

    int mpi_inited;
    MPI_Initialized(&mpi_inited);
    if (mpi_inited) {
	MPI_Comm_size(mpi_comm, &mpi_size);
	MPI_Comm_rank(mpi_comm, &mpi_rank);
    }

    hid_t fid, gid;
    hsize_t dims[2];
    H5T_class_t class_id;
    size_t type_size;
    char token[256];

    if (mpi_rank == 0) {
        fid = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
        if (fid < 0)  printf("Can not open %s\n", fn);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&fid, sizeof(hid_t), MPI_BYTE, 0, mpi_comm);
    }
    if (fid < 0) exit(1);

    if (mpi_rank == 0) {
	H5LTread_dataset_double(fid, "EWALD_L", ewald::L);
	H5LTread_dataset_double(fid, "TIME", &time);
	H5LTread_dataset_int(fid, "LT", &lt);

	// For compatibility
	if (H5LTfind_dataset(fid, "STRAIN")) {
	    H5LTread_dataset_double(fid, "STRAIN", &ewald::strain);
	}
	else if (H5LTfind_dataset(fid, "XSHEAR")) {
	    double xshear;
	    H5LTread_dataset_double(fid, "XSHEAR", &xshear);
	    ewald::strain = xshear/ewald::L[2];
	}

	H5LTread_dataset_double(fid, "SHRATE", &shRate);

	printf("time = %.3f\n", time);
	printf("lt = %d\n", lt);
	printf("L = %.3f %.3f %.3f\n", ewald::L[0], ewald::L[1], ewald::L[2]);
	printf("strain = %.3f\n", ewald::strain);
	printf("shear rate = %.3E\n", shRate);
    }
    if (mpi_size > 1) {
	MPI_Bcast(ewald::L, 3, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&time, 1, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&lt, 1, MPI_INT, 0, mpi_comm);
	MPI_Bcast(&ewald::strain, 1, MPI_DOUBLE, 0, mpi_comm);
	MPI_Bcast(&shRate, 1, MPI_DOUBLE, 0, mpi_comm);
    }

    int nvesicle =0;
    if (mpi_rank == 0) {
        if (H5LTfind_dataset(fid, "NVESICLE")) {
	    H5LTread_dataset_int(fid, "NVESICLE", &nvesicle);
        }
        printf("nvesicle = %d\n", nvesicle);
    }
    if (mpi_size > 1) {
	MPI_Bcast(&nvesicle, 1, MPI_INT, 0, mpi_comm);
    }

    // Read vesicles
    vesicles.resize(nvesicle);
    for (int ives = 0; ives < nvesicle; ives++) {
        Vesicle &vesicle = vesicles[ives];
        int nvert, nface;

	if (mpi_rank == 0) {
	    sprintf(token, "VESICLE%d", ives);
	    gid = H5Gopen(fid, token, H5P_DEFAULT);

	    H5LTget_dataset_info(gid, "X", dims, &class_id, &type_size);
	    nvert = dims[0];

	    H5LTget_dataset_info(gid, "F2V", dims, &class_id, &type_size);
	    nface = dims[0];
        }
	if (mpi_size > 1) {
	    MPI_Bcast(&nvert, 1, MPI_INT, 0, mpi_comm);
	    MPI_Bcast(&nface, 1, MPI_INT, 0, mpi_comm);
	}

        MArray<double,2> x(nvert,3);
	MArray<int,2> f2v(nface,3);

	if (mpi_rank == 0) {
	    H5LTread_dataset_double(gid, "X", x.data());
	    H5LTread_dataset_int(gid, "F2V", f2v.data());
        }
	if (mpi_size > 1) {
	    MPI_Bcast(x.data(), x.size(), MPI_DOUBLE, 0, mpi_comm);
	    MPI_Bcast(f2v.data(), f2v.size(), MPI_INT, 0, mpi_comm);
	}

	vesicle.verts.resize(nvert);
	vesicle.setCoords(x);

	vesicle.faces.resize(nface);
	vesicle.setConnectivities(f2v);

	if (mpi_rank == 0) H5Gclose(gid);
    } 

    if (mpi_rank == 0) H5Fclose(fid);
}


/* Write to file the centroid, velocity, and orientation of the vesicles
 * Arguments:
 *   fn -- file name
 */
void VeSys::writeVesicleCenter(const char *fn)
{
    int nves = numVesicles();
    if (nves <= 0) return;

    MArray<double,2> x(nves,3), v(nves,3);
    MArray<double,1> psi(nves);

    for (int ives = 0; ives < nves; ives++) {
        Vesicle &vesicle = vesicles[ives];

	m_dcopy(3, vesicle.center, &x(ives,0));
	vesicle.calcCentroidVel(vesicle.v, &v(ives,0));

	double psi_tmp, D_tmp;
	vesicle.shapeFactor(psi_tmp, D_tmp);
	psi(ives) = psi_tmp;
    }

    // Write to file
    hid_t fid, dset;
    int ierr;
    hsize_t dims[2];
    char token[256];

    // No error message
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    if (H5Fis_hdf5(fn) > 0) {
	fid = H5Fopen(fn, H5F_ACC_RDWR, H5P_DEFAULT);
    }
    else {
	fid = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    }

    if (! H5LTfind_dataset(fid, "NVESICLE")) {
        dims[0] = 1;
	H5LTmake_dataset_int(fid, "NVESICLE", 1, dims, &nves);
    }

    sprintf(token, "X%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nves;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, token, 2, dims, x.data());
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, x.data());
	H5Dclose(dset);
    }

    sprintf(token, "V%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nves;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, token, 2, dims, v.data());
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
	H5Dclose(dset);
    }
    
    sprintf(token, "PSI%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = nves;
	H5LTmake_dataset_double(fid, token, 1, dims, psi.data());
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, psi.data());
	H5Dclose(dset);
    }


    H5Fclose(fid);
}


/* Write the particle stress
 * Arguments:
 *   fn -- file name
 *   tau -- particle stress tensor */
void VeSys::writeStress(const char *fn, const double (*tau)[3])
{
    static bool first_write = true;
    FILE *file;

    if (first_write) {
        file = fopen(fn, "w");
        first_write = false;
    }
    else {
        file = fopen(fn, "a");
    }

    fprintf(file, "%12.5E %12.5E %12.5E %12.5E %12.5E\n",
                time, tau[0][2], tau[0][0], tau[1][1], tau[2][2]);

    fclose(file);
}


// Hematocrit, i.e. cell volume fraction
double VeSys::volFrac()
{
    double vol = 0.0;
    for (int ives = 0; ives < numVesicles(); ives++) {
        vol += vesicles[ives].vol;
    }

    return vol/ewald::vol;
}


// Maximum velocity perturbation
void VeSys::calcMaxVelPerturb(double *dvmax)
{
    m_dclear(3, dvmax);

    for (int ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];

	for (int ivert = 0; ivert < vesicle.numVerts(); ivert++) {
	    Point &vert = vesicle.verts[ivert];

	    double vbkg[3];
	    calcBkgVel(vert.x, vbkg);

	    FOR_D3 {
		double dv = vesicle.v(ivert,d) - vbkg[d];
		dvmax[d] = max(dvmax[d], fabs(dv));
	    }
	}
    } // ives
}


/* Solve the surface velocity and surface tensions
 * Note:
 *   The solution is stored in vesicle.v and vesicle.sigma */
void VeSys::solveVel()
{
    const int nvert_tot = vertList.size();

    //**********************************************************************
    // Prediction
    VecZeroEntries(rhs_vel);

    // Single-layer in rhs
    for (int p = 0, ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	if ( ! vesicle.isActive ) continue;

	int nvert = vesicle.numVerts();

	MArray<double,2> f(nvert,3);
	f = 0.0;
	f += vesicle.fbend;
	f += vesicle.ftens;

//	// debug
//	for (int ivert = 0; ivert < nvert; ivert++) {
//	    Point &vert = vesicle.verts[ivert];
//
//	    f(ivert,0) += 0.1*sin(vert.x[1]*2*M_PI*ewald::iL[1]);
//	    f(ivert,1) += 0.1*sin(vert.x[2]*2*M_PI*ewald::iL[2]);
//	    f(ivert,2) += 0.1*sin(vert.x[0]*2*M_PI*ewald::iL[0]);
//	}
//	// end debug

	// debug
	// push vesicles toward the middle of the domain in the z-direction
	double z = vesicle.center[2];
	const double zmin = 0.3*ewald::L[2];
	const double zmax = 0.7*ewald::L[2];

	z -= floor(z*ewald::iL[2])*ewald::L[2];
	if (z < zmin){
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        f(ivert,0) -= min(0.5, zmin - z);
	    }
	}
	else if (z > 0.7*ewald::L[2]) {
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        f(ivert,0) += min(0.5, z - zmax);
	    }
	}
	// end debug
	
	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = vesicle.verts[ivert];
	    m_dcopy(3, &f(ivert,0), vert.f);
	}
    }

    addBoundaryIntegral(-1.0/(8*M_PI), 0.0, rhs_vel);

    // Background velocity in rhs
    for (int ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	if ( ! vesicle.isPrivate ) continue;

	int nvert = vesicle.numVerts();

	MArray<double,2> v(nvert,3);
	v = 0.0;

	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = vesicle.verts[ivert];
	    calcBkgVel(vert.x, &v(ivert,0));
	}

	MArray<int,2> idx(nvert,3);
	for (int ivert = 0; ivert < nvert; ivert++) {
	    int p = 3*vesicle.verts[ivert].Gindx;
	    FOR_J3 idx(ivert,j) = p++;
	}

	VecSetValues(rhs_vel, 3*nvert, idx.data(), v.data(), ADD_VALUES);
    }

    VecAssemblyBegin(rhs_vel);
    VecAssemblyEnd(rhs_vel);

    // Initial guess
    VecZeroEntries(sol_vel);

    for (int ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	if ( ! vesicle.isPrivate ) continue;

	int nvert = vesicle.numVerts();

	MArray<int,2> idx(nvert,3);
	for (int ivert = 0; ivert < nvert; ivert++) {
	    int p = 3*vesicle.verts[ivert].Gindx;
	    FOR_J3 idx(ivert,j) = p++;
	}

	VecSetValues(sol_vel, 3*nvert, idx.data(), vesicle.v.data(), INSERT_VALUES);
    }

    VecAssemblyBegin(sol_vel);
    VecAssemblyEnd(sol_vel);

    // Solve BIE
    solveBie(Ts);

    // Copy solution
    MArray<double,2> v(nvert_tot,3);
    myVecScatter(scatter_vel, sol_vel, v.data(), INSERT_VALUES, SCATTER_FORWARD);

    for (int p = 0, ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	int nvert = vesicle.numVerts();

	m_dcopy(3*nvert, &v(p,0), vesicle.v.data());
	p += nvert;
    }

    //**********************************************************************
    // Projection
    const double tau = param::exist("AREA_RELAX_TIME") ? 
    		       param::getDoubleValue("AREA_RELAX_TIME") : 1.E99;

    // Target surface divergence
    MArray<double,2> v_old(nvert_tot,3), v_new(nvert_tot,3);
    MArray<double,1> divTar(nvert_tot);
    MArray<double,1> dlbd(nvert_tot);

    for (int p = 0, ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	int nvert = vesicle.numVerts();

	double s = (vesicle.areaTar/vesicle.area - 1.0)/tau;

	for (int ivert = 0; ivert < nvert; ivert++) {
	    FOR_D3 v_old(p,d) = vesicle.v(ivert,d);
	    divTar(p) = s*vesicle.vertArea(ivert);

	    p++;
	}
    } // ives

    dlbd = 0.0;
    projectVel(v_old, divTar, v_new, dlbd, true);

    // Update surface velocity and tension
    for (int p = 0, ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	int nvert = vesicle.numVerts();

	for (int ivert = 0; ivert < nvert; ivert++) {
	    FOR_D3 vesicle.v(ivert,d) = v_new(p,d);
	    vesicle.sigma(ivert) += dlbd(p);

	    p++;
	}
    }
}


void VeSys::timeInt()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    int nvert_tot = vertList.size();

    // Init
    for (int ives = 0; ives < numVesicles(); ives++) {
	Vesicle &vesicle = vesicles[ives];
	vesicle.v = 0.0;
	vesicle.sigma = 0.0;
    }

    Nt0 = lt;
    for (; lt <= Nt; lt++) {
	if (mpi_rank == 0) printf("lt = %9d time = %.5f\n", lt, time);
        double wtimeBgn = MPI_Wtime();

	for (int ives = 0; ives < numVesicles(); ives++) {
	    Vesicle &vesicle = vesicles[ives];

	    vesicle.updateGeometry();
	    vesicle.calcDoubleLayerJump();
	    vesicle.bendForce(vesicle.fbend);
	    vesicle.tensionForce(vesicle.sigma, vesicle.ftens);
	}

	domainDecomp();
	updateNeighborList();
	calcBieMat();
	solveVel();

	// Write the state before the surface evolves
	writeAll();

	// Diagnose
	double Ht = volFrac();

	double dvmax[3];
	calcMaxVelPerturb(dvmax);
	MPI_Bcast(dvmax, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double dA_min = FLT_MAX, dA_max = -FLT_MAX, dA_mean = 0.0, dA_cnt = 0.0;
	double div_min = FLT_MAX, div_max = -FLT_MAX, div_mean = 0.0, div_cnt = 0.0;

	for (int ives = 0; ives < numVesicles(); ives++) {
	    Vesicle &vesicle = vesicles[ives];
	    int nvert = vesicle.numVerts();

	    double dA = vesicle.area/vesicle.areaTar - 1.0; 
	    dA_min = min(dA_min, dA);
	    dA_max = max(dA_max, dA);
	    dA_mean += fabs(dA);
	    dA_cnt++;

	    MArray<double,1> div(nvert);
	    vesicle.velDiv(vesicle.v, div);
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        double tmp = div(ivert)/vesicle.vertArea(ivert);

		div_min = min(div_min, tmp);
		div_max = max(div_max, tmp);
		div_mean += fabs(tmp);
		div_cnt++;
	    }
	}

	dA_mean /= dA_cnt;
	div_mean /= div_cnt;

	if (mpi_rank == 0) {
	    printf("  vol frac   = %10.4f\n", Ht);
	    printf("  max perturb vel = %10.3f %10.3f %10.3f\n", dvmax[0], dvmax[1], dvmax[2]);
	    printf("  DlogA/Dt = %9.2E %9.2E %9.2E\n", div_min, div_max, div_mean);
	    printf("  area change = %9.2E %9.2E %9.2E\n", dA_min, dA_max, dA_mean);
	}

	if ( m_dmax(3, dvmax) > 10*max(1.0, fabs(shRate)*ewald::L[2])) {
	    if (mpi_rank == 0) {
	        printf("Warning: numerical instability detected\n");
	    }

	    break;
	}

	// Move the vesicles
	for (int ives = 0; ives < numVesicles(); ives++) {
	    Vesicle &vesicle = vesicles[ives];
	    int nvert = vesicle.numVerts();

	    MArray<double,2> vrelax(nvert,3);
	    vesicle.calcMeshRelaxVel(vrelax);

	    double UT[3], OMG[3];
	    vesicle.calcTransRotatVel(vesicle.v, UT, OMG);

	    for (int ivert = 0; ivert < nvert; ivert++) {
	        Point &vert = vesicle.verts[ivert];

		double v[3];
		FOR_D3 v[d] = vesicle.v(ivert,d) - UT[d];

		const double *nrml = &vesicle.vertNrml(ivert,0);
		double vn = m_ddot(3, nrml, v);

		FOR_D3 v[d] = UT[d] + vn*nrml[d] + vrelax(ivert,d);

		FOR_D3 vert.x[d] += Ts*v[d];
	    }

	    // Adjust volume
	    vesicle.updateAVC();
	    double Vtar = max(0.99*vesicle.vol, min(1.01*vesicle.vol, vesicle.volTar));
	    double s = pow(Vtar/vesicle.vol, 1.0/3.0);
	    for (int ivert = 0; ivert < nvert; ivert++) {
	        Point &vert = vesicle.verts[ivert];

		FOR_D3 {
		    double xx = vert.x[d] - vesicle.center[d];
		    vert.x[d] = vesicle.center[d] + s*xx;
	        }
	    }
	}

	ewald::strain += Ts*shRate;
	ewald::updateThreeAxes();

	// Prevent cell overlapping
	updateSourceGeometry();

	vector<NbrList*> nlists;
	nlists.push_back(&nlist_phys);
	collision::forceSeparation(vertList, nlists, 0.03);

	// Post-process
	rebox();
	syncCoord();
	time += Ts;

	// Time cost
        double wtimeEnd = MPI_Wtime();
	if (mpi_rank == 0) {
	    printf("    wall time = %7.2f s\n", wtimeEnd - wtimeBgn);
	}

	// Check whether to kill job
	int killjob = 0;
	if (mpi_rank == 0) {
	    const char fn[] = "killjob";
	    killjob = miscUtils::fileExists(fn);
	    if (killjob) {
	        printf("Found file %s, terminate job\n", fn);
	        remove(fn);
	    }
	}
	MPI_Bcast(&killjob, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (killjob) break;
    } // lt
}
