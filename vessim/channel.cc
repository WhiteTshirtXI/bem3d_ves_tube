#include "cxxheaders.h"
#include "mathfunc.h"
#include "ewald.h"
#include "stokesys.h"
#include "param.h"
#include "miscUtils.h"

void runjob(int argc, char **argv);

int main(int argc, char **argv)
{

    PetscInitialize(&argc, &argv, PETSC_NULL, PETSC_NULL);


    runjob(argc, argv);

    PetscFinalize();
}


void runjob(int argc, char **argv)
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    string tempstring, tempstring2;

    //cout << "Hello\n";

    // read parameters
    param::readFile("Input/channel.in");
    string fn_restart = param::getStringValue("RESTART");

    // read restart file
    StokeSys stks;

	//process command line overrides: for RESTART
	//Need to process RESTART override now before initialhook
	//I'm looking for command line information passed from the .pbs file that indicates if we want to use a different initial condition than what is listed in the input file.
	if(argc > 2)
	{
		tempstring = argv[1];
		if(tempstring.compare("-override") == 0) //0 means matches exactly
		{
			for(int ar = 2; ar < argc-1; ++ar)
			{
				tempstring = argv[ar];
				if(tempstring.compare("-RESTART") == 0)
				{
					fn_restart = argv[ar+1];
				}
			}
		}		
	}





    stks.readRestart(fn_restart.c_str());

    // Init Ewald sum
    if (param::exist("TIMING")) ewald::timing = true;
    ewald::alpha = square(0.2);
    ewald::tol = param::getDoubleValue("EWALD_TOL");
    ewald::initialHook();

    //cout << "Before initialhook\n";
    // Init the Stokes flow system
    stks.initialHook();
    //cout << "Afterinitialhook\n";

    // Set background velocity
    stks.vbkg[0] = param::getDoubleValue("VBKG", 0);
    stks.vbkg[1] = param::getDoubleValue("VBKG", 1);
    stks.vbkg[2] = param::getDoubleValue("VBKG", 2);

    // Set up cells
    double volRat = 1.0;  // need to specify the volume ratio in the channel.in
    if (param::exist("VOL_RAT")) volRat = param::getDoubleValue("VOL_RAT");
    double spontaneouscurv = 0;
    if (param::exist("H0")) spontaneouscurv = param::getDoubleValue("H0");



    for (int icell = 0; icell < stks.numCells(); ++icell) {
        Vesicle &cell = stks.cells[icell];

        // Material property
	//cell.ES = param::getDoubleValue("CELL_ES");
	//cell.ED = 100.0;
	cell.EB = param::getDoubleValue("CELL_EB");
	cell.H0 = spontaneouscurv;

	

	// Target volume
	cell.areaTar = 4*3.1415926535897932384626433832795; //Nondimensionalize length by radius of sphere of equivalent surface area                         
	cell.volTar = volRat * 4.0/3.0 * 3.1415926535897932384626433832795;
	
    }

    stks.viscRat = param::getDoubleValue("VISC_RAT");

    // Time integration
    stks.Ts = param::getDoubleValue("Ts");
    stks.Nt = param::getIntValue("Nt");
    //cout << mpi_rank << "Before initflowsolver\n";
    stks.initFlowSolver();
    //cout << mpi_rank << "After initflowsolver\n";


    //Command line overrides.  Optional: Override channel.in Input file directly from commands given in pbs script.  Or don't use additional command line flags after channel and ignore this.
    //Attempt to load output directory from file
    if (param::exist("OUTPUT_DIRECTORY")) stks.outfiledirectory = param::getStringValue("OUTPUT_DIRECTORY");
		else stks.outfiledirectory = "D";  //Default output directory

	if(argc > 2)
	{
		tempstring = argv[1];  //Check the first command line parameter.  If you see "-override" continue, other wise ignore everything else
		if(tempstring.compare("-override") == 0) //0 means matches exactly
		{
			for(int ar = 2; ar < argc-1; ar = ar+2)  //For each parameter after that
			{                                        //Look for an exact parameter name string match and update parameter accordingly
				tempstring = argv[ar];               
				if(tempstring.compare("-OUTPUT_DIRECTORY") == 0)  //Specify output directory from command line
				{
					stks.outfiledirectory = argv[ar+1];
				}
				else if(tempstring.compare("-VOL_RAT") == 0)  //Reduced volume ratio
				{
					tempstring2 = argv[ar+1];
					for(int icell = 0; icell < stks.numCells(); icell++)  //Change capillary number
					{					
						Vesicle &cell = stks.cells[icell];
						//cell.areaTar = 4*3.1415926535897932384626433832795; //Nondimensionalize length by radius of sphere of equivalent surface area                         
						cell.volTar = atof(tempstring2.c_str()) * 4.0/3.0 * 3.1415926535897932384626433832795;
					}
				}
				else if(tempstring.compare("-H0") == 0)  //Spontaneous curvature
				{
					tempstring2 = argv[ar+1];
					for(int icell = 0; icell < stks.numCells(); icell++)  //Change capillary number
					{					
						Vesicle &cell = stks.cells[icell];
						cell.H0 = atof(tempstring2.c_str());
					}
				}
				else if(tempstring.compare("-ca") == 0)  //Capillary number.  reciprocal is multiplier for EB for special cell.  
				{
					tempstring2 = argv[ar+1];
					
					for(int icell = 0; icell < stks.numCells(); icell++)  //Change capillary number
					{					
						Vesicle &cell = stks.cells[icell];
						cell.EB = cell.EB*1./atof(tempstring2.c_str());
						//cell.ED = cell.ED*1./atof(tempstring2.c_str());
					}
				}
				else if(tempstring.compare("-VISC_RAT") == 0)  //Capillary number.  reciprocal is multiplier for ES and EB for special cell.  
				{
					tempstring2 = argv[ar+1];
					stks.viscRat = atof(tempstring2.c_str());
				}				
			}
		}
	}
	//cout << "Hello again\n";

    miscUtils::mkdir(stks.outfiledirectory.c_str());
    //cout << "Entering timeINT\n";
 
   // The real work of this function occurs in stokesys.cc's timeInt(), which in turn references stokesys_flowsolver.cc's solveFlow().
   stks.timeInt();

    // Finalize
    ewald::finalHook();
}
