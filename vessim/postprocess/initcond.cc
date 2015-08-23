#include "cxxheaders.h"
#include "ewald.h"
#include "stokesys.h"
#include "mathfunc.h"
#include "param.h"
#include "miscUtils.h"

StokeSys stks;
Cell cellRef;

int main()
{
    MPI_Init(NULL, NULL);

    const int ranseed = 12345;
    srand(ranseed);
    cout << "ranseed = " << ranseed << " -> " << rand01() << endl;

    // Wall
    stks.walls.resize(1);
    Wall *wall = &stks.walls[0];
    wall->readHDF5("Input/channel_16x9.h5");

    // Periodic domain
    {
	double xmin[3], xmax[3];
	wall->getCoordRange(xmin, xmax);

	using ewald::L;
	L[0] = xmax[0] - xmin[0];
	L[1] = xmax[1] - xmin[1];
	L[2] = 30./2.82;		// height of the channel

	ewald::vol = L[0]*L[1]*L[2];
	ewald::iVol = 1.0/ewald::vol;
	printf("Domain size = %.3f %.3f %.3f\n", L[0], L[1], L[2]); 

	for (int ivert = 0; ivert < wall->numVerts(); ivert++) {
	    Point &vert = wall->verts[ivert];
	    FOR_I3 vert.x[i] -= xmin[i];
	}
    }

    // Reference cell
    cellRef.makeBiConcave(8);
    cellRef.setInternalPointers();
    cellRef.updateGeometry();
    double szCell[3];
    {
	double xmin[3], xmax[3];
	cellRef.getCoordRange(xmin, xmax);
	FOR_I3 szCell[i] = xmax[i] - xmin[i];
    }


    // Maximum cell occupation
    int N0, N1, N2;
    double szbox[3];
    {
	N0 = (int)floor(ewald::L[0]/(1.1*szCell[0]));
	N1 = (int)floor(ewald::L[1]/(1.05*szCell[1]));
	N2 = (int)floor(ewald::L[2]/(1.1*szCell[2]));

	szbox[0] = ewald::L[0]/N0;
	szbox[1] = ewald::L[1]/N1;
	szbox[2] = ewald::L[2]/N2;

        printf("nbox = %d %d %d\n", N0, N1, N2);
    }

    double Htmax = N0*N1*N2*cellRef.vol/ewald::vol;
    printf("Max Ht = %.3f\n", Htmax);

    // Place cells
    double Ht_tar;
    printf("Ht =? ");
    cin >> Ht_tar;
    Ht_tar = max(0.01, min(Ht_tar, Htmax));
    int ncell_tar = (int)round(ewald::vol*Ht_tar/cellRef.vol);

    MArray<bool,3> occupied(N0,N1,N2);
    occupied = false;

    while (stks.numCells() < ncell_tar) {
        int i0, i1, i2;
	double xc[3];

	i0 = (int)floor(rand01()*N0);
	i1 = (int)floor(rand01()*N1);
	i2 = (int)floor(rand01()*N2);
	if (occupied(i0,i1,i2)) continue;

	// Create a new cell at the fcc lattice sites
	xc[0] = (i0 + 0.5)*szbox[0] + 0.01*(rand01() - 1.0)*szCell[0];
	xc[1] = (i1 + 0.5)*szbox[1] + 0.01*(rand01() - 1.0)*szCell[1];
	xc[2] = (i2 + 0.5)*szbox[2] + 0.01*(rand01() - 1.0)*szCell[2];

	stks.cells.push_back(Cell());

	Cell *cell = &stks.cells.back();
	cell->copyMeshFrom(cellRef);

	cell->cellRef = new Cell();
	cell->cellRef->copyMeshFrom(cellRef);

	// Scale and translate
	for (int ivert = 0; ivert < cell->numVerts(); ivert++) {
	    Point &vert = cell->verts[ivert];
	    FOR_I3 vert.x[i] += xc[i];
	}

	// mark the block as occupied
	occupied(i0,i1,i2) = true;
    } // while


    printf("Number of cells placed = %d\n", stks.numCells());
    printf("Ht = %.3f\n", stks.numCells()*cellRef.vol/ewald::vol);

    // Write results
    stks.time = 0.0;
    stks.lt = 0;
    FOR_I3 stks.vbkg[i] = 0.0;

    miscUtils::mkdir("D");

    char fn[1024];
    sprintf(fn, "%s%6.6d%s", "D/restart", stks.lt, ".dat");
    stks.writeRestart(fn);

    sprintf(fn, "%s%6.6d%s", "D/cell", stks.lt, ".dat");
    stks.writeCells(fn);

    sprintf(fn, "%s%6.6d%s", "D/wall", stks.lt, ".dat");
    stks.writeWalls(fn);

    // Finalize
    MPI_Finalize();
}
