// Read a restart file and then move platelets to center so that z0 < z < z1
#include "cxxheaders.h"
#include "stokesys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "param.h"
#include "miscUtils.h"


int main(int argc, char **argv)
{
    // Init
    MPI_Init(&argc, &argv);

    // Init random number
    const int ranseed = 12345;
    srand(ranseed);
    cout << "ranseed = " << ranseed << " -> " << rand01() << endl;

    // Init system
    StokeSys stks;
    stks.readRestart(argv[1]);
//    stks.readRestart("D/restart020000.dat");
    stks.initialHook();

    for (int icell = 0; icell < stks.numCells(); icell++) {
        Cell &cell = stks.cells[icell];
	cell.updateGeometry();
    }
    stks.rigids.resize(0);

    // Read reference platelet and rotate till the shortest axis
    // align with z-axis
    Rigid rigidRef;
    //rigidRef.readHDF5("Input/platelet.h5");
    rigidRef.readDATsingle("D/seansphere.dat");
    rigidRef.setInternalPointers();
    rigidRef.updateGeometry();

    // Rotate the platelet so that its axis of rotation is along x-axis
    { 
	double xmax[3] = { -FLT_MAX };
	double xmin[3] = { FLT_MAX };
	double laxis[3];

	for (int ivert = 0; ivert < rigidRef.numVerts(); ivert++) {
	    Point &vert = rigidRef.verts[ivert];

	    FOR_I3 {
	        double tmp = m_ddot(3, rigidRef.paxis[i], vert.x);
	        xmax[i] = max(xmax[i], tmp);
	        xmin[i] = min(xmin[i], tmp);
	    }
	}

	FOR_I3 laxis[i] = xmax[i] - xmin[i];

	int ix = 0;
	for (int i = 1; i < 3; i++) {
	    if (laxis[i] < laxis[ix]) ix = i;
	}

	// Rotate the rigid
	double A[3][3];
	double ex[3] = {1.0, 0.0, 0.0};
	calcRotateMatrix(rigidRef.paxis[ix], ex, A);

	for (int ivert = 0; ivert < rigidRef.numVerts(); ivert++) {
	    Point &vert = rigidRef.verts[ivert];

	    double xx[3];
	    FOR_I3 xx[i] = vert.x[i] - rigidRef.center[i];

	    FOR_I3 {
	        vert.x[i] = rigidRef.center[i];
		FOR_J3 vert.x[i] += A[i][j]*xx[j];
	    }
	}

	rigidRef.updateGeometry();
    } // end 


    // Determine the size of the box
    double szbox[3];
    int N0, N1, N2;

    {
        using ewald::L;

	double xmin[3], xmax[3];

	rigidRef.getCoordRange(xmin, xmax);

	szbox[0] = 1.0*(xmax[0] - xmin[0]);
	szbox[1] = 1.2*(xmax[1] - xmin[1]);
	szbox[2] = 1.2*(xmax[2] - xmin[2]);

	N0 = (int)floor(L[0]/szbox[0]);		
	N1 = (int)floor(L[1]/szbox[1]);		
	N2 = (int)floor(L[2]/szbox[2]);
	//if (N2%2 == 0) N2--;
	//N2 = 10; //Manual override that  may be useful

	szbox[0] = L[0]/N0;
	szbox[1] = L[1]/N1;	
	szbox[2] = L[2]/N2;

	cout << "xwidth = " << xmax[0]-xmin[0] << endl;
	cout << "ywidth = " << xmax[1]-xmin[1] << endl;
	cout << "zwidth = " << xmax[2]-xmin[2] << endl;
	cout << "Horizontal slices = " << N2 << endl;
    }


    // Find all boxes that are occupied by red cells
    MArray<bool,3> occupied(N0,N1,N2);
    occupied = false;

    for (int icell = 0; icell < stks.numCells(); icell++) {
        Cell &cell = stks.cells[icell];

	double xmin[3], xmax[3];
	int imin[3], imax[3];

	cell.getCoordRange(xmin, xmax);

	FOR_I3 {
	    imin[i] = (int)floor(xmin[i]/szbox[i]);
	    imax[i] = (int)ceil(xmax[i]/szbox[i]);
	}

	MArray<bool,3> mask(IndexRange(imin[0], imax[0]), 
			    IndexRange(imin[1], imax[1]),
			    IndexRange(imin[2], imax[2]));
	mask = false;

	// Make all the boxes that collide with the cell surface
	for (int iface = 0; iface < cell.numFaces(); iface++) {
	    double xtri[3][3];
	    FOR_I3
	    FOR_J3 {
	        xtri[i][j] = cell.faces[iface].vert[i]->x[j];
	    }

	    double xtrimin[3], xtrimax[3];
	    FOR_I3 {
	        xtrimin[i] = min(xtri[0][i], min(xtri[1][i], xtri[2][i]));
	        xtrimax[i] = max(xtri[0][i], max(xtri[1][i], xtri[2][i]));
	    }

	    for (int i0 = (int)floor(xtrimin[0]/szbox[0]); i0 < (int)ceil(xtrimax[0]/szbox[0]); i0++)
	    for (int i1 = (int)floor(xtrimin[1]/szbox[1]); i1 < (int)ceil(xtrimax[1]/szbox[1]); i1++)
	    for (int i2 = (int)floor(xtrimin[2]/szbox[2]); i2 < (int)ceil(xtrimax[2]/szbox[2]); i2++) {
		double LB[3] = { szbox[0]*i0, szbox[1]*i1, szbox[2]*i2 };
		double UB[3] = { szbox[0]*(i0+1), szbox[1]*(i1+1), szbox[2]*(i2+1) };

		if ( tri_box_collide(xtri, LB, UB) ) {
		    mask(i0,i1,i2) = true;
		    // debug
		    // printf("%3d: %3d %3d %3d\n", icell, i0, i1, i2);
		    // end debug
	        }
	    }
	}


	// Use a flood method to fill the 3D mask
	for (int i1 = imin[1]; i1 < imax[1]; i1++)
	for (int i2 = imin[2]; i2 < imax[2]; i2++) {
	    int i0bgn, i0end;

	    for (i0bgn = imin[0]; i0bgn < imax[0]; i0bgn++) {
	        if (mask(i0bgn,i1,i2)) break;
	    }

	    for (i0end = imax[0]; i0end > imin[0]; i0end--) {
	        if (mask(i0end-1,i1,i2)) break;
	    }

	    for (int i0 = i0bgn; i0 < i0end; i0++) {
	        int j0 = modulo(i0, N0);
		int j1 = modulo(i1, N1);
		int j2 = i2;
		occupied(j0,j1,j2) = true;
	    }
	}
    }



    // Further reduce the availability of boxes
    for (int i0 = 0; i0 < N0; i0++)
    for (int i1 = 0; i1 < N1; i1++)
    for (int i2 = 0; i2 < N2; i2++) {
	// Only add rigid at the center line
	//if ( i2 == 0 || i2 == N2-1) occupied(i0,i1,i2) = true;

	// Do not place platelets in two boxes next to each other in x-direction
	if (! occupied(i0,i1,i2) ) {
	    occupied((i0+1)%N0, i1, i2) = true;
	}
    }

    // Find all available boxes
    vector<int> i0_avail, i1_avail, i2_avail;
    for (int i0 = 0; i0 < N0; i0++)
    for (int i1 = 0; i1 < N1; i1++)
    for (int i2 = 0; i2 < N2; i2++) {
        if (!occupied(i0,i1,i2)) {
	    i0_avail.push_back(i0);
	    i1_avail.push_back(i1);
	    i2_avail.push_back(i2);
	}
    }

    int nrigid_avail = i0_avail.size();
    printf("Number of available rigid sites = %d\n", nrigid_avail);

    int nrigid;
    printf("Number of platelets =? ");
    cin >> nrigid;
    nrigid = min(nrigid, nrigid_avail);

    stks.rigids.resize(nrigid);
    for (int irigid = 0; irigid < nrigid; irigid++) {
	Rigid &rigid = stks.rigids[irigid];
	rigid.copyMeshFrom(rigidRef);

	// Pick up an available box
	int p = (int)floor(rand01()*nrigid_avail);
	int i0 = i0_avail[p];
	int i1 = i1_avail[p];
	int i2 = i2_avail[p];

	// Remove the available box
	i0_avail.erase(i0_avail.begin() + p);
	i1_avail.erase(i1_avail.begin() + p);
	i2_avail.erase(i2_avail.begin() + p);
	nrigid_avail--;


	double xx[3];
	xx[0] = (i0 + 0.5)*szbox[0] - rigidRef.center[0];
	xx[1] = (i1 + 0.5)*szbox[1] - rigidRef.center[1];
	xx[2] = (i2 + 0.5)*szbox[2] - rigidRef.center[2];

	for (int ivert = 0; ivert < rigid.numVerts(); ivert++) {
	    Point &vert = rigid.verts[ivert];
	    FOR_I3 vert.x[i] += xx[i];
	}
    } // irigid

    // Remove probes and tracers
    stks.nprb = 0;
    stks.ntrac = 0;

    // Write restart and tecplot file
    // Reset time step counter
    stks.lt = 0;
    stks.time = 0.0;

    miscUtils::mkdir("D");

    char fn[1024];
    sprintf(fn, "%s%6.6d%s", "D/restart", stks.lt, ".dat");
    stks.writeRestart(fn);
    printf("New restart file: %s\n", fn);

    sprintf(fn, "%s%6.6d%s", "D/cell", stks.lt, ".dat");
    stks.writeCells(fn);
    printf("New cell file: %s\n", fn);

    sprintf(fn, "%s%6.6d%s", "D/rigid", stks.lt, ".dat");
    stks.writeRigids(fn);
    printf("New platelet file: %s\n", fn);

    sprintf(fn, "%s%6.6d%s", "D/wall", stks.lt, ".dat");
    stks.writeWalls(fn);
    printf("New wall file: %s\n", fn);

    // Finalize
    MPI_Finalize();
}
