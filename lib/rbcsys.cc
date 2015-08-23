#include "cxxheaders.h"
#include "rbcsys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "quadrature.h"
#include "param.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "miscUtils.h"
#include "debugfunc.h"
#include "collision.h"


/* RbcSys::RbcSys 
 * RbcSys::~RbcSys 
 * RbcSys::initialHook 
 * RbcSys::updateGeometry 
 * RbcSys::updateSourceList
 * RbcSys::updateTargetList 
 * RbcSys::updateNeighborList
 * RbcSys::timeInt
 * RbcSys::rebox
 * RbcSys::syncCoord
 * RbcSys::calcStress
 * RbcSys::writeAll 
 * RbcSys::writeCells 
 * RbcSys::writeRestart 
 * RbcSys::readRestart 
 * RbcSys::bcastRestart
 * RbcSys::writeProbe
 * RbcSys::writeCellCenter
 * RbcSys::writeStress 
 * RbcSys::cellVolFraction
 * RbcSys::calcMaxVelPerturb
 * RbcSys::calcAreaChange */


RbcSys::RbcSys()
{
    matSL = NULL;
    matDL = NULL;
}


RbcSys::~RbcSys()
{
}

/* -- Index meshes, vertices, and faces
 * -- Build vertex lists and face lists
 * -- Compute geometries of reference cells */
void RbcSys::initialHook()
{
    // Build mesh list and index meshes
    for (int i = 0; i < numCells(); i++) {
        RedCell &cell = cells[i];
        meshList.push_back(&cell);

	cell.Gindx = i;
	cell.setInternalPointers();
	cell.noPeriodicBoundaries();
        cell.buildEdgeList();
        cell.buildControlPoints();

	RedCell *cellRef = cell.cellRef;
	if (cellRef) {
	    cellRef->setInternalPointers();
	    cellRef->noPeriodicBoundaries();
	    cellRef->buildEdgeList();
	    cellRef->updateGeometry();
	}
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

    //======================================================================
    // Allocate f and v array for each mesh
    for (int icell = 0; icell < numCells(); icell++) {
        RedCell &cell = cells[icell];

	int nvert = cell.numVerts();
	cell.f.resize(nvert,3);
	cell.v.resize(nvert,3);
    }
}



void RbcSys::updateGeometry()
{
    for (int i = 0; i < numCells(); i++) {
        cells[i].updateGeometry();
    }
}


void RbcSys::updateSourceList()
{
    // Clear mesh active flag
    vector<Mesh*> &meshes = meshList;
    for (int i = 0; i < meshes.size(); i++) {
        meshes[i]->active = false;
    }

    // Physical
    if (ewald::phys::active) {
        SrcList &slist = slist_phys;
        slist.build(faceList);

	for (int i = 0; i < slist.faces.size(); i++) {
	    Mesh *mesh = slist.faces[i]->mesh;
	    if (mesh) mesh->active = true;
	}
    }

    // Fourier
    if (ewald::four::active) {
        vector<Tri*> &slist = slist_four;
	ewald::four::setSources(faceList, slist);

	for (int i = 0; i < slist.size(); i++) {
	    Mesh *mesh = slist[i]->mesh;
	    if (mesh) mesh->active = true;
	}
    }
}


void RbcSys::updateTargetList()
{
    if (ewald::four::active) {
        ewald::four::setTargets(vertList, tlist_four);
    }
}


void RbcSys::updateNeighborList()
{
    if (ewald::phys::active) {
	NbrList &nlist = nlist_phys;
	nlist.build(vertList, slist_phys);
    }
}


// Compute the single- and double-layer matrices
void RbcSys::calcBieMat()
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


// Time integration
void RbcSys::timeInt()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    Nt0 = lt;

    for (; lt <= Nt; lt++) {

	if (mpi_rank == 0) printf("lt = %9d time = %.5f\n", lt, time);
        double wtimeBgn = MPI_Wtime();

	// Domain decomposition
	updateGeometry();
	updateSourceList();
	updateTargetList();
	updateNeighborList();

        // Compute cell surface force
	for (int icell = 0; icell < numCells(); icell++) {
	    RedCell &cell = cells[icell];

	    cell.area_E2 = cell.area_E1;
	    cell.area_E1 = cell.area_E0;
	    cell.area_E0 = cell.area/cell.areaTar - 1.0;
	    cell.sigma += cell.KP*(cell.area_E0 - cell.area_E1) 
			+ cell.KI* cell.area_E0
			+ cell.KD*(cell.area_E0 - 2*cell.area_E1 + cell.area_E2);

	    // f = elastic + bend + viscous + tension
	    cell.f = 0.0;

	    int nvert = cell.numVerts();
	    MArray<double,2> ftmp(nvert,3);
	    cell.elasticForce(ftmp);
	    cell.f += ftmp;

	    cell.bendForce(ftmp);
	    cell.f += ftmp;

	    cell.viscousForce(cell.v, ftmp);
	    cell.f += ftmp;

	    cell.tensionForce(cell.sigma, ftmp);
	    cell.f += ftmp;

//	    // debug
//	    // Stir the system
//	    for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
//	        Point &vert = cell.verts[ivert];
//		cell.f(ivert,0) += sin(2*M_PI*ewald::iL[1]*vert.x[1]);
//		cell.f(ivert,1) += sin(2*M_PI*ewald::iL[2]*vert.x[2]);
//		cell.f(ivert,2) += sin(2*M_PI*ewald::iL[0]*vert.x[0]);
//	    }
//	    // end debug
	}

	calcBieMat();
	calcResidual_Bie(rhs_vel);
	solveBie(Ts);

	// Copy solution
	int nvert_tot = vertList.size();
	MArray<double,2> v(nvert_tot,3);
	myVecScatter(scatter_vel, sol_vel, v.data(), INSERT_VALUES, SCATTER_FORWARD);

	for (int p = 0, icell = 0; icell < numCells(); icell++) {
	    RedCell &cell = cells[icell];
	    int nvert = cell.numVerts();

	    m_dcopy(3*nvert, &v(p,0), cell.v.data());
	    p += nvert;
	}

	
	// Output the state before system evolves
        if (lt != Nt0) writeAll();

	// Diagnose
	double Ht = cellVolFraction();

	double dvmax[3];
	calcMaxVelPerturb(dvmax);
	MPI_Bcast(dvmax, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	double dAmin, dAmax, dAmean;
	calcAreaChange(dAmin, dAmax, dAmean);

	if (mpi_rank == 0) {
	    if(numCells() > 0) {
		printf("    cell vol frac   = %10.3f\n", Ht);
		printf("    max perturb vel = %10.3f %10.3f %10.3f\n", 
			    dvmax[0], dvmax[1], dvmax[2]);
		printf("    area change (min,max,mean) = %10.3f %10.3f %10.3f\n", 
			    dAmin, dAmax, dAmean);
            }
	}

	if ( m_dmax(3, dvmax) > 100*max(1.0, fabs(shRate)*ewald::L[2])) {
	    if (mpi_rank == 0) {
	        printf("Numerical instability detected, stop\n");
	    }

	    break;
	}

	// Evolve the system
	for (int icell = 0; icell < numCells(); icell++) {
	    RedCell &cell = cells[icell];
	    for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	        m_daxpy(3, Ts, &cell.v(ivert,0), cell.verts[ivert].x);
	    }

	    // Adjust volume
	    cell.calcAreaAndVolume(cell.area, cell.vol);
	    double Vtar = max(0.99*cell.vol, 
				min(1.01*cell.vol, cell.volTar));
	    double s = pow(Vtar/cell.vol, 1.0/3.0);
	    for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	        Point &vert = cell.verts[ivert];

		FOR_D3 {
		    double xx = vert.x[d] - cell.center[d];
		    vert.x[d] = cell.center[d] + s*xx;
	        }
	    }
	    cell.updateGeometry();
	}

	ewald::xshear += Ts*shRate*ewald::L[2];
	ewald::xshear -= nearbyint(ewald::xshear*ewald::iL[0])*ewald::L[0];

	// Prevent cell overlapping
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


// Rebox cells and rigid particles
void RbcSys::rebox()
{
    using ewald::L;
    using ewald::iL;

    vector<Mesh*> &meshes = meshList;

    for (int imesh = 0; imesh < meshes.size(); imesh++) {
        Mesh *mesh = meshes[imesh];

        // Determine translating distance
        double xmin[3], xmax[3], xc[3];

        mesh->getCoordRange(xmin, xmax);
	FOR_D3 xc[d] = 0.5*(xmin[d] + xmax[d]);

	bool move_mesh = false;
	FOR_D3 {
            double eps = 0.01*L[d];
	    if (xc[d] < -eps || xc[d] > L[d] + eps) move_mesh = true;
	}

	if (move_mesh) {
	    double origin[3];
	    FOR_D3 origin[d] = 0.5*ewald::L[d];

	    double xcnew[3];
	    FOR_D3 xcnew[d] = xc[d];
	    ewald::to_CloseBy(origin, xcnew);

	    double xx[3];
	    FOR_D3 xx[d] = xcnew[d] - xc[d];

	    for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
		Point &vert = mesh->verts[ivert];
		FOR_I3 vert.x[i] += xx[i];
	    }
	}
    } // imesh
}


// Synchronize cell and rigid particle coordinates
void RbcSys::syncCoord()
{
    // Assemble global coordinate arrays
    vector<Mesh*> &meshes = meshList;

    int npoint = 0;
    for (int imesh = 0; imesh < meshes.size(); imesh++)
        npoint += meshes[imesh]->numVerts();

    double (*x)[3] = new double[npoint][3];

    for (int imesh = 0, p = 0; imesh < meshes.size(); imesh++) {
        Mesh *mesh = meshes[imesh];
	for (int ivert = 0; ivert < mesh->numVerts(); ivert++) {
	    m_dcopy(3, mesh->verts[ivert].x, x[p]);
	    p++;
	}
    }

    // Broadcast from root node
    MPI_Bcast(*x, 3*npoint, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Update vertex coordinates
    for (int imesh = 0, p = 0; imesh < meshes.size(); imesh++) {
        Mesh *mesh = meshes[imesh];
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
 *   tau -- added stress, normalized by (cell vol) * (shear rate)  */
void RbcSys::calcStress(double (*tau)[3])
{
    // Init
    m_dclear(9, *tau);

    for (int icell = 0; icell < numCells(); icell++) {
        RedCell &cell = cells[icell];

	for (int ifa = 0; ifa < cell.faces.size(); ifa++) {
	    Tri &face = cell.faces[ifa];

	    double xtri[3][3], ftri[3][3], vtri[3][3]; 
	    for (int l = 0; l < 3; l++) {
		int ivert = face.ivert[l];
	        m_dcopy(3, cell.verts[ivert].x, xtri[l]);
		m_dcopy(3, &cell.f(ivert,0), ftri[l]);
		m_dcopy(3, &cell.v(ivert,0), vtri[l]);
	    }

	    Quad2D &Q = quadrature::select_rule_2d("TRI_3");
	    for (int iq = 0; iq < Q.n(); iq++) {
	        double s = Q.x(iq), t = Q.y(iq);
		double r = 1.0 - s - t;

	        double xq[3], fq[3], vq[3], dA;
		FOR_I3 {
		    xq[i] = r*xtri[0][i] + s*xtri[1][i] + t*xtri[2][i] - cell.center[i];
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
    } // icell

    // normalize by (total cell volume)*(shear rate)
    double cvol = 0.0;
    for (int icell = 0; icell < numCells(); icell++) {
        RedCell &cell = cells[icell];
        cvol += cell.vol;
    }
    m_dscal(9, 1.0/(cvol*shRate), *tau);
} 


// Write everything
void RbcSys::writeAll()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    const char FN_FMT[] = "%s%6.6d%s";
    char token[256], fn[256];
    int nout;
    const int LT_CHUNK = 10000;

    // Cells
    strcpy(token, "CELL_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/cell", lt, ".dat");
	if (mpi_rank == 0) writeCells(fn);
    }


    // Restart
    strcpy(token, "RESTART_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
	sprintf(fn, FN_FMT, "D/restart", lt, ".dat");
	if (mpi_rank == 0) writeRestart(fn);
    }

    // Center of cells
    strcpy(token, "CELL_CENTER_OUT");
    nout = param::exist(token)? param::getIntValue(token) : -1;
    if (nout > 0 && lt%nout == 0) {
        sprintf(fn, FN_FMT, "D/cell_center", (lt/LT_CHUNK)*LT_CHUNK, ".dat");
        if (mpi_rank == 0) writeCellCenter(fn);
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


/* Write all cells to file
 * Arguments:
 *   fn -- file name */
void RbcSys::writeCells(const char *fn)
{
    if (numCells() <= 0) return;

    FILE *file = fopen(fn, "w");
//    fprintf(file, "variables = x, y, z, u, v, w\n");
    fprintf(file, "variables = x, y, z\n");

    for (int icell = 0; icell < numCells(); icell++) {
        RedCell &cell = cells[icell];
	int nvert = cell.numVerts();
	int nface = cell.numFaces();
	fprintf(file, "zone N=%d E=%d F=FEPOINT ET=TRIANGLE\n", nvert, nface);
        
	// coordinates
	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = cell.verts[ivert];
	    fprintf(file, "%10.3e %10.3e %10.3e\n", vert.x[0], vert.x[1], vert.x[2]);
        }

	// connectivity
	for (int iface = 0; iface < nface; iface++) {
	    Tri &face = cell.faces[iface];
	    fprintf(file, "%d  %d  %d\n", face.ivert[0]+1, face.ivert[1]+1, face.ivert[2]+1);
        }
    }

    fclose(file);
}



void RbcSys::writeRestart(const char *fn)
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
    H5LTmake_dataset_double(fid, "XSHEAR", 1, dims, &ewald::xshear);

    dims[0] = 1;
    H5LTmake_dataset_double(fid, "SHRATE", 1, dims, &shRate);

    int ncell = numCells();
    if (ncell > 0) {
	dims[0] = 1;
	H5LTmake_dataset_int(fid, "NCELL", 1, dims, &ncell);
    }

    // Cells 
    for (int icell = 0; icell < numCells(); icell++) {
        RedCell &cell = cells[icell];
	int nvert = cell.numVerts();
	int nface = cell.numFaces();

	MArray<double,2> x(nvert,3);
	MArray<int,2> f2v(nface,3);
	cell.getCoords(x);
	cell.getConnectivities(f2v);

        sprintf(token, "CELL%d", icell);
	hid_t gid = H5Gcreate(fid, token, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	dims[0] = nvert;
	dims[1] = 3;
        H5LTmake_dataset_double(gid, "X", 2, dims, x.data());

	if (cell.cellRef != NULL) {
	    cell.cellRef->getCoords(x);
	    dims[0] = nvert;
	    dims[1] = 3;
	    H5LTmake_dataset_double(gid, "XREF", 2, dims, x.data());
	}

	dims[0] = nface;
	dims[1] = 3;
	H5LTmake_dataset_int(gid, "F2V", 2, dims, f2v.data());

	// Surface tension
	dims[0] = 1;
	H5LTmake_dataset_double(gid, "CELL_SIGMA", 1, dims, &cell.sigma);
	H5LTmake_dataset_double(gid, "CELL_AREA_E0", 1, dims, &cell.area_E0);
	H5LTmake_dataset_double(gid, "CELL_AREA_E1", 1, dims, &cell.area_E1);
	H5LTmake_dataset_double(gid, "CELL_AREA_E2", 1, dims, &cell.area_E2);

	H5Gclose(gid);
    }

    H5Fclose(fid);
}


/* Read restart file in HDF5 format
 * Arguments:
 *  fn -- file name */
void RbcSys::readRestart(const char *fn)
{
    const MPI_Comm comm = MPI_COMM_WORLD;
    int rank;
    MPI_Comm_rank(comm, &rank);

    hid_t fid, gid;
    hsize_t dims[2];
    H5T_class_t class_id;
    size_t type_size;
    char token[256];

    if (rank == 0) {
	fid = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (fid < 0)  printf("Can not open %s\n", fn);
    }
    MPI_Bcast(&fid, sizeof(hid_t), MPI_BYTE, 0, comm);
    if (fid < 0) exit(1);

    if (rank == 0) {
	H5LTread_dataset_double(fid, "EWALD_L", ewald::L);
	H5LTread_dataset_double(fid, "TIME", &time);
	H5LTread_dataset_int(fid, "LT", &lt);
	H5LTread_dataset_double(fid, "XSHEAR", &ewald::xshear);
	H5LTread_dataset_double(fid, "SHRATE", &shRate);

	printf("time = %.3f\n", time);
	printf("lt = %d\n", lt);
	printf("L = %.3f %.3f %.3f\n", ewald::L[0], ewald::L[1], ewald::L[2]);
	printf("xshear = %.3f\n", ewald::xshear);
	printf("shear rate = %.3E\n", shRate);
    }
    MPI_Bcast(ewald::L, 3, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&time, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&lt, 1, MPI_INT, 0, comm);
    MPI_Bcast(&ewald::xshear, 1, MPI_DOUBLE, 0, comm);
    MPI_Bcast(&shRate, 1, MPI_DOUBLE, 0, comm);

    int ncell =0;
    if (rank == 0) {
	if (H5LTfind_dataset(fid, "NCELL")) {
	    H5LTread_dataset_int(fid, "NCELL", &ncell);
        }
        printf("ncell = %d\n", ncell);
    }
    MPI_Bcast(&ncell, 1, MPI_INT, 0, comm);


    // Read cells
    cells.resize(ncell);

    for (int icell = 0; icell < ncell; icell++) {
        RedCell &cell = cells[icell];
        int nvert, nface, xref_flag;

	if (rank == 0) {
	    sprintf(token, "CELL%d", icell);
	    gid = H5Gopen(fid, token, H5P_DEFAULT);

	    H5LTget_dataset_info(gid, "X", dims, &class_id, &type_size);
	    nvert = dims[0];

	    H5LTget_dataset_info(gid, "F2V", dims, &class_id, &type_size);
	    nface = dims[0];

	    xref_flag = H5LTfind_dataset(gid, "XREF");
	}
	MPI_Bcast(&nvert, 1, MPI_INT, 0, comm);
	MPI_Bcast(&nface, 1, MPI_INT, 0, comm);
	MPI_Bcast(&xref_flag, 1, MPI_INT, 0, comm);

        MArray<double,2> x(nvert,3), xref(nvert,3);
	MArray<int,2> f2v(nface,3);

	if (rank == 0) {
	    H5LTread_dataset_double(gid, "X", x.data());
	    H5LTread_dataset_int(gid, "F2V", f2v.data());
	    if (xref_flag) H5LTread_dataset_double(gid, "XREF", xref.data());
        }
	MPI_Bcast(x.data(), x.size(), MPI_DOUBLE, 0, comm);
	MPI_Bcast(f2v.data(), f2v.size(), MPI_INT, 0, comm);
	if (xref_flag) {
	    MPI_Bcast(xref.data(), xref.size(), MPI_DOUBLE, 0, comm);
	}

	cell.verts.resize(nvert);
	cell.setCoords(x);

	cell.faces.resize(nface);
	cell.setConnectivities(f2v);

        if (xref_flag) {
	    cell.cellRef = new RedCell();

	    cell.cellRef->verts.resize(nvert);
	    cell.cellRef->setCoords(xref);

	    cell.cellRef->faces.resize(nface);
	    cell.cellRef->setConnectivities(f2v);
	}

	// Surface tension
	if (rank == 0) {
	    if (H5LTfind_dataset(gid, "CELL_SIGMA")) {
		H5LTread_dataset_double(gid, "CELL_SIGMA", &cell.sigma);
	    }

	    if (H5LTfind_dataset(gid, "CELL_AREA_E0")) {
		H5LTread_dataset_double(gid, "CELL_AREA_E0", &cell.area_E0);
		H5LTread_dataset_double(gid, "CELL_AREA_E1", &cell.area_E1);
		H5LTread_dataset_double(gid, "CELL_AREA_E2", &cell.area_E2);
	    }
	}
	MPI_Bcast(&cell.sigma, 1, MPI_DOUBLE, 0, comm);
	MPI_Bcast(&cell.area_E0, 1, MPI_DOUBLE, 0, comm);
	MPI_Bcast(&cell.area_E1, 1, MPI_DOUBLE, 0, comm);
	MPI_Bcast(&cell.area_E2, 1, MPI_DOUBLE, 0, comm);

	if (rank == 0) H5Gclose(gid);
    } 

    if (rank == 0) H5Fclose(fid);
}


void RbcSys::writeCellCenter(const char *fn)
{
    int ncell = numCells();
    if (ncell <= 0) return;

    double (*x)[3], (*v)[3];
    x = new double[ncell][3];
    v = new double[ncell][3];

    for (int icell = 0; icell < ncell; icell++) {
        RedCell &cell = cells[icell];

	FOR_I3 x[icell][i] = cell.center[i];
	cell.calcCentroidVel(cell.v, v[icell]);
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

    if (! H5LTfind_dataset(fid, "NCELL")) {
        dims[0] = 1;
	H5LTmake_dataset_int(fid, "NCELL", 1, dims, &ncell);
    }

    sprintf(token, "X%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = ncell;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, token, 2, dims, *x);
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *x);
	H5Dclose(dset);
    }

    sprintf(token, "V%6.6d", lt);
    if (! H5LTfind_dataset(fid, token)) {
        dims[0] = ncell;
	dims[1] = 3;
	H5LTmake_dataset_double(fid, token, 2, dims, *v);
    }
    else {
        dset = H5Dopen(fid, token, H5P_DEFAULT);
	H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, *v);
	H5Dclose(dset);
    }
    
    H5Fclose(fid);
    delete [] x;
    delete [] v;
}


/* Write the add stress */
void RbcSys::writeStress(const char *fn, const double (*tau)[3])
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
double RbcSys::cellVolFraction()
{
    double vol = 0.0;
    for (int icell = 0; icell < numCells(); icell++) {
        vol += cells[icell].vol;
    }

    return vol/ewald::vol;
}


// Maximum velocity perturbation
void RbcSys::calcMaxVelPerturb(double *dvmax)
{
    m_dclear(3, dvmax);

    for (int icell = 0; icell < numCells(); icell++) {
	RedCell &cell = cells[icell];

	for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	    Point &vert = cell.verts[ivert];

	    double vbkg[3];
	    calcBkgVel(vert.x, vbkg);

	    FOR_I3 {
		double dv = cell.v(ivert,i) - vbkg[i];
		dvmax[i] = max(dvmax[i], fabs(dv));
	    }
	}
    } // icell
}


// Max area change
void RbcSys::calcAreaChange(double &dAmin, double &dAmax, double &dAmean)
{
    // Init
    dAmin = 1.E10;
    dAmax = -1.E10;

    double dAsum = 0.0, Asum = 0.0;

    for (int icell = 0; icell < numCells(); icell++) {
	RedCell &cell = cells[icell];
	RedCell *cellRef = cell.cellRef;
	if (cellRef == NULL) continue;

	double dA = cell.area/cellRef->area - 1.0;
	dAmin = min(dAmin, dA);
	dAmax = max(dAmax, dA);

	dAsum += fabs(dA);
	Asum += cellRef->area;
    } // icell

    dAmean = dAsum/Asum;
    if (dAmin > dAmax) dAmin = dAmax = 0.0;
}
