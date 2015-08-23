// Calculate the pressure gradient
#include "cxxheaders.h"
#include "stokesys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "param.h"

string in_dir("D15/ht10/ca1/");
string out_dir("D15/anawall/");

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    param::readFile("Input/channel.in");
    param::readFile("Input/analyze_wall.in");

    string fn_restart = "D15/ht10/ca1/restart969000.dat";//param::getStringValue("RESTART");
    int lt_min = 968900;//param::getIntValue("LT_MIN");
    int lt_max = 983200;//param::getIntValue("LT_MAX");
    int lt_step = 10000;//param::getIntValue("LT_STEP");

    if (param::exist("IN_DIR")) in_dir = "D15/ht10/ca1/"; //param::getStringValue("IN_DIR") + "/";
    if (param::exist("OUT_DIR")) out_dir = "D15/anawall/"; //param::getStringValue("OUT_DIR") + "/";

    // Read restart
    StokeSys stks;
    stks.readRestart(fn_restart.c_str());
    stks.initialHook();
    for (int iwall = 0; iwall < stks.numWalls(); iwall++) 
    {
		Wall &wall = stks.walls[iwall];
		wall.updateGeometry();
    }

    ewald::alpha = 1.0;
    ewald::tol = 0.99;
    ewald::set_prms();

    // Read wall files
    Wall &wall = stks.walls[0];

    int nrecMax = (lt_max - lt_min)/lt_step + 1;
    int nrec = 0;
    MArray<double,2> pgrad(nrecMax, 3);

    for (int lt = lt_min; lt <= lt_max; lt += lt_step) {
	char fn[256];
	const int BUFF_SIZE = 1024;
	char buff[BUFF_SIZE];

	sprintf(fn, "%s%s%6.6d%s", in_dir.c_str(), "wall", lt, ".dat");
	FILE *file = fopen(fn, "r");
	if (! file) continue;

	fgets(buff, BUFF_SIZE, file);
	fgets(buff, BUFF_SIZE, file);
	for (int ivert = 0; ivert < wall.numVerts(); ivert++) {
	    double xtmp[3], ftmp[3];
	    fgets(buff, BUFF_SIZE, file);
	    sscanf(buff, "%lf %lf %lf %lf %lf %lf", 
	    		&xtmp[0], &xtmp[1], &xtmp[2], 
	    		&ftmp[0], &ftmp[1], &ftmp[2]); 

	    FOR_I3 wall.verts[ivert].f[i] = ftmp[i];
	    
	    cout << "Vert " << ivert << " x: " << xtmp[0] << ' ' << xtmp[1] << ' ' << xtmp[2] << endl;
	    cout << "Vert " << ivert << " f: " << ftmp[0] << ' ' << ftmp[1] << ' ' << ftmp[2] << endl;

	}
	fclose(file);
	printf("Read %s\n", fn);

	stks.calcPressGrad(&pgrad(nrec,0));
	cout << "Pressure gradient for nrec = " << nrec << ": " << pgrad(nrec,0) << ' ' << pgrad(nrec,1) << ' ' << pgrad(nrec,2) << endl;
	nrec++;
    }

    double pgrad_MEAN[3] = {0.0};
    double pgrad_RMS[3] = {0.0};

    for (int irec = 0; irec < nrec; irec++) {
        FOR_J3 {
	    pgrad_MEAN[j] += pgrad(irec,j);
	    pgrad_RMS[j] += square(pgrad(irec,j));
        }
    }

    FOR_I3 pgrad_MEAN[i] /= nrec;
    FOR_I3 pgrad_RMS[i] = sqrt(pgrad_RMS[i]/nrec - square(pgrad_MEAN[i]));

    printf("Number of records = %d\n", nrec);
    printf("Mean pgrad = %13.3E %13.3E %13.3E\n", 
    		pgrad_MEAN[0], pgrad_MEAN[1], pgrad_MEAN[2]);
    printf("Wall shear rate = %13.3E\n", pgrad_MEAN[0]*0.5*ewald::L[2]);
    printf("RMS pgrad  = %13.3E %13.3E %13.3E\n", 
    		pgrad_RMS[0], pgrad_RMS[1], pgrad_RMS[2]);

    if (param::exist("VEL_PROFILE")) {
        char fn_out[512];
	sprintf(fn_out, "%s%s", out_dir.c_str(), "Poiseuille.dat");
	printf("Write velocity profile to %s\n", fn_out);

	FILE *fout = fopen(fn_out, "w");

	int n = 128;
	fprintf(fout, "# dp/dx = %13.3E\n", pgrad_MEAN[0]);
	fprintf(fout, "VARIABLES = z, u\n");
	fprintf(fout, "ZONE I=%d F=POINT T=\"%s\"\n", n+1, "Poiseuille");
	for (int i = 0; i <= n; i++) {
	    double z = i*ewald::L[2]/n;
	    double u = 0.5*pgrad_MEAN[0]*z*(ewald::L[2] - z);
	    fprintf(fout, "%15.5E %15.5E\n", z, u);
	}

	fclose(fout);
    }

    MPI_Finalize();
}
