#include "cxxheaders.h"
#include "ewald.h"
#include "veswall.h"
#include "mathfunc.h"
#include "param.h"
#include "miscUtils.h"

VesWall stks;
Vesicle cellRef;

int main()
{
    MPI_Init(NULL, NULL);

	//Commented out: No randomness in this version
    //const int ranseed = 12345;
    //srand(ranseed);
    //cout << "ranseed = " << ranseed << " -> " << rand01() << endl;

    // Load a Wall
    stks.walls.resize(1);
    Wall *wall = &stks.walls[0];
    
    //HDF5 reader:
    //wall->readHDF5("Input/cylinder005.h5");

    //Read from DAT file:
    wall->readDATsingle("Input/wall20_128.dat");


    // Periodic domain
    {
	double xmin[3], xmax[3];
	wall->getCoordRange(xmin, xmax);

	//Domain dimensions
	using ewald::L;
	L[0] = xmax[0] - xmin[0]; 
	L[1] = 3;  // Tube: make this bigger than the actual channel (e.g. add +0.5 or multiply by 1.05, etc.)
	//L[2] = xmax[2] - xmin[2];  
	L[2] = 3;		// height of the channel.  Must be slightly greater than diameter of wall mesh

	ewald::vol = L[0]*L[1]*L[2];
	ewald::iVol = 1.0/ewald::vol;
	printf("Domain size = %.3f %.3f %.3f\n", L[0], L[1], L[2]); 

	for (int ivert = 0; ivert < wall->numVerts(); ivert++) {
	    Point &vert = wall->verts[ivert];
	    FOR_I3 vert.x[i] -= xmin[i];
	}
    }

	double xcen, ycen, zcen; // this is hardcoded by Andrew, but you would want to replace this with cin, cout etc. to specify the cell's center.
        xcen = 10;   //Center of 5-long channel
		ycen = 1.; //Center of 2.2-wide channel
		zcen = 1.; //Center of 2.2-high channel


//    // Reference cell
//    cellRef.makeBiConcave(8);
//    cellRef.setInternalPointers();
//    cellRef.updateGeometry();
//    double szCell[3];
//    {
//	double xmin[3], xmax[3];
//	cellRef.getCoordRange(xmin, xmax);
//	FOR_I3 szCell[i] = xmax[i] - xmin[i];
//    }



	    //Now add the cell.
		//For this version it's a rigid, but it will be a cell soon
		Vesicle newves;
                //newves.readDATsingle("/work/03201/spann/blood95/85/conf83/v1eb1/D/vesicle000000.dat");
		newves.readDATsingle("Input/ves65eqlow.dat");  //Set the mesh you want to load here (pregenerated meshes for vesicles if you want to use them)
		// reads a tecplot-compliant file rather than a hdf5-compliant file
		// files should be centered at 0,0
		newves.setInternalPointers();
		newves.updateGeometry();
		for (int ivert = 0; ivert < newves.numVerts(); ivert++) 
		{
			Point &vert = newves.verts[ivert];
			vert.x[0] += xcen;
			vert.x[1] += ycen;
			vert.x[2] += zcen;
		}
		
		
		newves.updateGeometry();
		stks.vesicles.push_back(newves);







    printf("Number of vesicles in system = %d\n", stks.numVesicles());
    //printf("Ht = %.3f\n", stks.numCells()*cellRef.vol/ewald::vol);

    // Write results
    stks.time = 0.0;
    stks.lt = 0;
    FOR_I3 stks.vbkg[i] = 0.0;

    miscUtils::mkdir("D");  //Output directory to write files to

    char fn[1024];
    sprintf(fn, "%s%6.6d%s", "D/restart", stks.lt, ".dat");  //Retrieve your file at D/restart000000.dat
    stks.writeRestart(fn);

    sprintf(fn, "%s%6.6d%s", "D/vesicle", stks.lt, ".dat");
    stks.writeVesiclesplain(fn);

    sprintf(fn, "%s%6.6d%s", "D/wall", stks.lt, ".dat");
    stks.writeWalls(fn);

    // Finalize
    MPI_Finalize();
}

