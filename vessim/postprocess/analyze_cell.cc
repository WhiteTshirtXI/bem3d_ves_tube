#include "cxxheaders.h"
#include "ewald.h"
#include "stokesys.h"
#include "param.h"
#include "mathfunc.h"
#include "miscUtils.h"

StokeSys stks;

double Ts;
int lt_min, lt_max, lt_step;
string in_dir("D15sean/ht10/ca1/");
string out_dir("D15sean/anacell/");
const int CHRLEN = 1024;

// Slice in z
const int nz = 24;
const double &Lz = ewald::L[2];

int nrec, nrigid;
// 1st index -- record number
// 2nd index -- cell number
// 3rd index -- coordinates
MArray<double,3> xcell, vcell;
MArray<int,1> steps;

// Volume fraction and stresslet, and residual force
MArray<double,1> vf_mean(nz); 
MArray<double,1> Sxz_mean(nz), Sxx_mean(nz), Syy_mean(nz), Szz_mean(nz);
MArray<double,1> fx_mean(nz), fy_mean(nz), fz_mean(nz);

// Volume fraction and stresslet, local definition
MArray<double,1> VolSum(nz); 
MArray<double,1> SxzSum(nz), SxxSum(nz), SyySum(nz), SzzSum(nz);

// Inclination angle
MArray<int,1> npsi(nz);
MArray<double,1> psi_mean(nz), psi_rms(nz);

double incline_angle(const Mesh &mesh);

int main(int argc, char **argv)
{
    param::readFile("Input/channel.in");
    //param::readFile("Input/analyze_cell.in");
    string tempstring;

    Ts = param::getDoubleValue("Ts");
    lt_min = 3000;//param::getIntValue("LT_MIN");
    lt_max = 30000;//param::getIntValue("LT_MAX");
    lt_step = 100;//param::getIntValue("LT_STEP");
    string restart_file = "D15sean/ht10/ca05/restart000000.dat";//param::getStringValue("RESTART");
    //tempstring = argv[0];
    //in_dir = tempstring;
    if (param::exist("IN_DIR")) in_dir = param::getStringValue("IN_DIR") + "/";
    if (param::exist("OUT_DIR")) out_dir = param::getStringValue("OUT_DIR") + "/";

    bool traj_flag = true; //param::exist("TRAJECTORY");
    bool volf_flag = param::exist("VOL_FRAC");
    bool stress_flag = param::exist("STRESS");
    bool psi_flag = param::exist("INCLINE_ANGLE");

    stks.readRestart(restart_file.c_str());

    // Only need to set up Ewald parameters, not the whole ewald::initialHook
    // just need some positive number here
    ewald::alpha = 1.0;		
    ewald::tol = 0.99;
    ewald::set_prms();

    for (int icell = 0; icell < stks.numCells(); ++icell) {
        Cell &cell = stks.cells[icell];

        // Material property
        cell.ES = param::getDoubleValue("CELL_ES");
        cell.ED = 100.0;
        cell.EB = 3.3E-3*cell.ES;
    }
    stks.initialHook();

    int ncell = stks.numCells();
    int nrecMax = (lt_max - lt_min)/lt_step + 1;
    xcell.resize(nrecMax,ncell,3);
    vcell.resize(nrecMax,ncell,3);
    steps.resize(nrecMax);

    // Init
    nrec = 0;

    vf_mean = 0.0;

    Sxz_mean = 0.0;
    Sxx_mean = 0.0;
    Syy_mean = 0.0;
    Szz_mean = 0.0;

    fx_mean = 0.0;
    fy_mean = 0.0;
    fz_mean = 0.0;

    SxzSum = 0.0;
    SxxSum = 0.0;
    SyySum = 0.0;
    SzzSum = 0.0;

    VolSum = 0.0;

    npsi = 0;
    psi_mean = 0.0;
    psi_rms = 0.0;

    for (int lt = lt_min; lt <= lt_max; lt += lt_step) {
	char fn_in[CHRLEN], line[CHRLEN];
	FILE *fin, *fout;

        sprintf(fn_in, "%s%s%6.6d%s", in_dir.c_str(), "cell", lt, ".dat");
	fin = fopen(fn_in, "r");
	if (! fin) continue;
	printf("Read %s\n", fn_in);

	// Skip the file head
	fgets(line, CHRLEN, fin);

	for (int icell = 0; icell < stks.cells.size(); ++icell) {
	    Cell &cell = stks.cells[icell];

	    fgets(line, CHRLEN, fin);	// skip zone header

	    // Read coordinates
	    for (int ivert = 0; ivert < cell.numVerts(); ivert++) {
	        fgets(line, CHRLEN, fin);	
	        Point &vert = cell.verts[ivert];
	        sscanf(line, "%lf %lf %lf", &vert.x[0], &vert.x[1], &vert.x[2]);
	    }
	    for (int iface = 0; iface < cell.numFaces(); iface++) {
	        fgets(line, CHRLEN, fin);	
	    }
	}
	fclose(fin);

	// Compute geometry and surface force
	for (int icell = 0; icell < stks.cells.size(); icell++) {
	    Cell &cell = stks.cells[icell];
	    cell.updateGeometry();
	    cell.calcSurfForce(cell.f);
	} 

	// Trajectory
	for (int icell = 0; icell < ncell; icell++) {
	    Cell &cell = stks.cells[icell];
	    FOR_I3 xcell(nrec,icell,i) = cell.center[i];
	}

        // Volume fraction
	if (volf_flag) {
	    MArray<double,1> vf(nz);
	    stks.calcCellVolFrac(0.0, Lz, nz, vf.data());
	    vf_mean += vf;
        }

	// Particle stress
	if (stress_flag) {
	    MArray<double,1> Sxz(nz), Sxx(nz), Syy(nz), Szz(nz), Vol(nz);
	    MArray<double,1> fx(nz), fy(nz), fz(nz);

	    stks.calcStress(0.0, Lz, nz, 
	    		Sxz.data(), Sxx.data(), Syy.data(), Szz.data(),
			fx.data(), fy.data(), fz.data());

	    Sxz_mean += Sxz;
	    Sxx_mean += Sxx;
	    Syy_mean += Syy;
	    Szz_mean += Szz;

	    fx_mean += fx;
	    fy_mean += fy;
	    fz_mean += fz;

	    stks.calcStressLocal(0.0, Lz, nz, 
	    		Sxz.data(), Sxx.data(), Syy.data(), Szz.data(), 
			Vol.data());

	    SxzSum += Sxz;
	    SxxSum += Sxx;
	    SyySum += Syy;
	    SzzSum += Szz;
	    VolSum += Vol;
	}

	// Inclination angle
	if (psi_flag) {
	    for (int icell = 0; icell < ncell; icell++) {
		Cell &cell = stks.cells[icell];

		double psi = incline_angle(cell);
		int iz = (int)floor(cell.center[2]/Lz*nz);

		npsi(iz)++;
		psi_mean(iz) += psi;
		psi_rms(iz) += psi*psi;
	    } 
	}

	steps(nrec) = lt;
	nrec++;	
    } //lt

    printf("Number of files = %d\n", nrec);

    // Normalization
    vf_mean *= 1.0/nrec;

    Sxz_mean *= 1.0/nrec;
    Sxx_mean *= 1.0/nrec;
    Syy_mean *= 1.0/nrec;
    Szz_mean *= 1.0/nrec;

    fx_mean *= 1.0/nrec;
    fy_mean *= 1.0/nrec;
    fz_mean *= 1.0/nrec;

    for (int iz = 0; iz < nz; iz++) {
        double in = 1.0/(npsi(iz) + 1.E-10);
        psi_mean(iz) *= in;
	psi_rms(iz) = sqrt(psi_rms(iz)*in - psi_mean(iz)*psi_mean(iz));
    }


    // Start writing output file 
    miscUtils::mkdir(out_dir.c_str());

    // Cell trajectory
    if (traj_flag) {
        char fn_out[CHRLEN];
	FILE *fout;

	sprintf(fn_out, "%s%s%6.6d%s%6.6d%s",
		    out_dir.c_str(), "cell_traj.", 
		    lt_min, "_to_", lt_max, ".dat");

	fout = fopen(fn_out, "w");
	fprintf(fout, "variables = t, x, y, z\n");
	for (int icell = 0; icell < ncell; icell++) {
	    fprintf(fout, "zone\n");
	    for (int irec = 0; irec < nrec; irec++) {
	        fprintf(fout, "%15.5E %15.5E %15.5E %15.5E\n", 
			steps(irec)*Ts,
			xcell(irec,icell,0),
			xcell(irec,icell,1),
			xcell(irec,icell,2));
	    }
	}
	fclose(fout);

	printf("Cell trajectory written to %s\n", fn_out);
    }


    // Cell velocity
    if (traj_flag) {
        char fn_out[CHRLEN];
	FILE *fout;

	double *u_mean = new double[nz];
	std::fill_n(u_mean, nz, 0.0);

	int *cnt = new int[nz];
	std::fill_n(cnt, nz, 0);


	double umeanall = 0;
	int reccountall = 0;
	for (int irec = 0; irec < nrec-1; irec++)
	for (int icell = 0; icell < ncell; icell++) {
	    double dx = xcell(irec+1, icell, 0) 
	              - xcell(irec, icell, 0);
	    dx -= nearbyint(dx*ewald::iL[0])*ewald::L[0];

	    double z = 0.5*(xcell(irec+1,icell,2) + xcell(irec,icell,2));
	    int iz = (int)(nz*z/Lz);

	    u_mean[iz] += dx/(Ts*lt_step);
	    umeanall += dx/(Ts*lt_step);
	    reccountall++;
	    cnt[iz]++;
	}

	for (int iz = 0; iz < nz; iz++) {
	    if (cnt[iz] > 0) u_mean[iz] /= cnt[iz];
	}


	sprintf(fn_out, "%s%s%6.6d%s%6.6d%s",
		out_dir.c_str(), "cell_vel.", lt_min, "_to_", lt_max, ".dat");
	fout = fopen(fn_out, "w");
	fprintf(fout, "variables = z, u\n");
	fprintf(fout, "zone");
	for (int iz = 0; iz < nz; iz++) {
	    double z = (iz+0.5)*Lz/nz;
	    fprintf(fout,  "%15.5e %15.5e\n", z, u_mean[iz]);
	}
	fclose(fout);

	printf("Cell x-vel written to %s\n", fn_out);

	umeanall = umeanall / double(reccountall);
	cout << "Umeanall = " << umeanall << endl;
	
    }

    if (volf_flag) {
        char fn_out[CHRLEN];
	sprintf(fn_out, "%s%s%6.6d%s%6.6d%s",
		    out_dir.c_str(), "cell_Ht.", 
		    lt_min, "_to_", lt_max, ".dat");
			    
	FILE *fout = fopen(fn_out, "w");
	fprintf(fout, "variables = z, vf\n");
	fprintf(fout, "zone\n");
	for (int i = 0; i < nz; i++) {
	    double z = (i+0.5)*Lz/nz;
	    fprintf(fout, "%15.5e %15.5e \n", z, vf_mean(i));
	}
	fclose(fout);

	printf("Cell vol fraction written to %s\n", fn_out);
    }

    if (stress_flag) {
	// Global stresslet
        char fn_out[CHRLEN];
	sprintf(fn_out, "%s%s%6.6d%s%6.6d%s",
		    out_dir.c_str(), "stress.", 
		    lt_min, "_to_", lt_max, ".dat");
			    
	FILE *fout = fopen(fn_out, "w");
	fprintf(fout, "variables = z, Sxz, Sxx, Syy, Szz, fx, fy, fz\n");
	fprintf(fout, "zone\n");
	for (int i = 0; i < nz; i++) {
	    double z = (i+0.5)*Lz/nz;
	    fprintf(fout, " %12.5e", z);
	    fprintf(fout, " %12.5e %12.5e %12.5e %12.5e", 
	    		Sxz_mean(i), Sxx_mean(i), Syy_mean(i), Szz_mean(i));
	    fprintf(fout, " %12.5e %12.5e %12.5e", 
	    		fx_mean(i), fy_mean(i), fz_mean(i));
	    fprintf(fout, "\n");
	}
	fclose(fout);
	printf("Particle stress written to %s\n", fn_out);

	// Local stresslet
	double Vslab = ewald::vol/nz;
	sprintf(fn_out, "%s%s%6.6d%s%6.6d%s",
		    out_dir.c_str(), "localStress.", 
		    lt_min, "_to_", lt_max, ".dat");
	fout = fopen(fn_out, "w");
	fprintf(fout, "variables = z, Sxz, Sxx, Syy, Szz, Sxz_per_cell, Sxx_per_Cell, Syy_per_cell, Szz_per_Cell, Ht\n");
	fprintf(fout, "zone\n");
	for (int i = 0; i < nz; i++) {
	    double z = (i+0.5)*Lz/nz;
	    fprintf(fout, " %12.5e", z);
	    fprintf(fout, " %12.5e %12.5e %12.5e %12.5e", 
	    		SxzSum(i)/(nrec*Vslab), 
	    		SxxSum(i)/(nrec*Vslab), 
	    		SyySum(i)/(nrec*Vslab), 
	    		SzzSum(i)/(nrec*Vslab) );

	    fprintf(fout, " %12.5e %12.5e %12.5e %12.5e", 
	    		SxzSum(i)/(VolSum(i) + 1.E-10), 
	    		SxxSum(i)/(VolSum(i) + 1.E-10), 
	    		SyySum(i)/(VolSum(i) + 1.E-10), 
	    		SzzSum(i)/(VolSum(i) + 1.E-10) );

	    fprintf(fout, " %12.5e", VolSum(i)/(nrec*Vslab));

	    fprintf(fout, "\n");
	}
	fclose(fout);
	printf("Local particle stress written to %s\n", fn_out);
    }


    if (psi_flag) {
        char fn_out[512];
	sprintf(fn_out, "%s%s%6.6d%s%6.6d%s",
		    out_dir.c_str(), "cell_psi.", lt_min, "_to_", lt_max, ".dat");
			    
	FILE *fout = fopen(fn_out, "w");
	fprintf(fout, "variables = z, psiMean, psiRms\n");
	fprintf(fout, "zone\n");
	for (int i = 0; i < nz; i++) {
	    double z = (i+0.5)*Lz/nz;
	    fprintf(fout, " %12.5e %12.5e %12.5e\n", 
	    		z, psi_mean(i)/M_PI, psi_rms(i)/M_PI );
	}
	fclose(fout);

	printf("Cell incline angle written to %s\n", fn_out);
    }
}


/* Calculate the deviation of axis of symmetry away from the 
 * z-axis
 * Argument
 *   mesh -- */
double incline_angle(const Mesh &mesh)
{
    // Length of the vesicle along the principle axes
    double Lmin[3], Lmax[3], L[3];
    FOR_D3 {
        Lmin[d] = FLT_MAX;
        Lmax[d] = -FLT_MAX;
    }

    for (int ivert = 0; ivert < mesh.numVerts(); ivert++) {
        const Point &vert = mesh.verts[ivert];
        FOR_D3 {
            double Ltmp = m_ddot(3, vert.x, mesh.paxis[d]);
            Lmin[d] = min(Lmin[d], Ltmp);
            Lmax[d] = max(Lmax[d], Ltmp);
        }
    }

    FOR_D3 L[d] = Lmax[d] - Lmin[d];

    // Identify the axis that is parallel to y-axis
    int iy = 0;
    for (int i = 1; i < 3; i++) {
        if (fabs(mesh.paxis[i][1]) > fabs(mesh.paxis[iy][1])) iy = i;
    }

    // For the rest two axes
    // ix is the longer axis and iz is the shorter one
    int ix = (iy+1)%3;
    for (int i = 0; i < 3; i++) {
        if (i == iy) continue;
        if (L[i] > L[ix]) ix = i;
    }

    int iz = 3 - ix - iy;

    return acos(fabs(mesh.paxis[iz][2]));
}
