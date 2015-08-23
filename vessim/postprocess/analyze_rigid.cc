#include "cxxheaders.h"
#include "stokesys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "param.h"
#include "miscUtils.h"
#include "hdf5.h"
#include "hdf5_hl.h"

StokeSys stks;

const double EPS = 1.E-10;

string in_dir("D/");
string out_dir("tmp/");
int lt_min, lt_max, lt_step;
double Ts;

int nrec, nrigid;
MArray<double,3> xrigid, vrigid;	// (irec, irigid, dim)
MArray<int,1> step;

void readRigid();
void trajectory();
void density();
void velProf();
void correlation();
void margination();

int main(int argc, char **argv)
{
    // Parameters
    param::readFile("Input/channel.in");
    param::readFile("Input/analyze_rigid.in");

    Ts = param::getDoubleValue("Ts");
    lt_min = param::getIntValue("LT_MIN");
    lt_max = param::getIntValue("LT_MAX");
    lt_step = param::getIntValue("LT_STEP");
    if (param::exist("IN_DIR")) in_dir = param::getStringValue("IN_DIR") + "/";
    if (param::exist("OUT_DIR")) out_dir = param::getStringValue("OUT_DIR") + "/";

    string restart_file = param::getStringValue("RESTART");
    string tempstring;  //Process command line overrides
    string tempstring2;
    
  	if(argc > 2)    //Specify output directory from command line, overriding analyze_rigid.in text file input
	{
		tempstring = argv[1];
		if(tempstring.compare("-override") == 0) //0 means matches exactly
		{
			for(int ar = 2; ar < argc-1; ++ar)
			{
				tempstring = argv[ar];
				if(tempstring.compare("-RESTART") == 0)
				{
					restart_file = argv[ar+1];
				}
				else if(tempstring.compare("-IN_DIR") == 0)
				{
					in_dir = argv[ar+1];
				}
				else if(tempstring.compare("-OUT_DIR") == 0)  
				{
					out_dir = argv[ar+1];
				}
				else if(tempstring.compare("-LT_MIN") == 0)  
				{
					tempstring2 = argv[ar+1];
					lt_min = atoi(tempstring2.c_str());
				}
				else if(tempstring.compare("-LT_MAX") == 0)  
				{
					tempstring2 = argv[ar+1];
					lt_max = atoi(tempstring2.c_str());
				}
				else if(tempstring.compare("-LT_STEP") == 0) 
				{
					tempstring2 = argv[ar+1];
					lt_step = atoi(tempstring2.c_str());
				}
				else if(tempstring.compare("-TS") == 0)  
				{
					tempstring2 = argv[ar+1];
					Ts = atof(tempstring2.c_str());
				}		
			}
		}
	}


    
    
    
    
    
    
    stks.readRestart(restart_file.c_str());
    // Only need to set up Ewald parameters, not the whole ewald::initialHook
    ewald::alpha = 1.0;	// just need some positive number here
    ewald::tol = 0.99;
    ewald::set_prms();

    // Alloc memory
    int nrecMax = (lt_max - lt_min)/lt_step + 1;
    nrigid = stks.numRigids();
    xrigid.resize(nrecMax,nrigid,3);
    vrigid.resize(nrecMax,nrigid,3);
    step.resize(nrecMax);

    readRigid();

    // Write output file 
    miscUtils::mkdir(out_dir.c_str());
    if (param::exist("TRAJECTORY")) trajectory();
    if (param::exist("DENSITY")) density();
    if (param::exist("VEL_PROFILE")) velProf();
    if (param::exist("CORRELATION")) correlation();
    if (param::exist("MARGINATION")) margination();
}


void readRigid()
{
    const int LT_CHUNK = 10000;
    char fn[256], token[256];
    hid_t fid = -1;
    double *xbuf = new double[3*nrigid];
    double *vbuf = new double[3*nrigid];

    // Turn off error message
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    nrec = 0;
    for (int lt = lt_min; lt <= lt_max; lt += lt_step) {
	static int lt_file_old = -1;
	int lt_file = (lt/LT_CHUNK)*LT_CHUNK;

	// Open a new file when necessary
	if (lt_file > lt_file_old) {
	    if (fid > 0) H5Fclose(fid);

	    sprintf(fn, "%s%s%6.6d%s", in_dir.c_str(), "rigid_center", lt_file, ".dat");
	    printf("    %s\n", fn);

	    fid = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
	    lt_file_old = lt_file;
	}
	if (fid <= 0) continue;

	sprintf(token, "%s%6.6d", "X", lt);
	if (H5LTfind_dataset(fid, token)) {
	    H5LTread_dataset_double(fid, token, xbuf);
	} else {
	    continue;
	}

	sprintf(token, "%s%6.6d", "V", lt);
	if (H5LTfind_dataset(fid, token)) {
	    H5LTread_dataset_double(fid, token, vbuf);
	} else {
	    continue;
	}

	step(nrec) = lt;
	m_dcopy(3*nrigid, xbuf, &xrigid(nrec,0,0));
	m_dcopy(3*nrigid, vbuf, &vrigid(nrec,0,0));
	nrec++;
    } // lt

    printf("Number of records = %d\n", nrec);

    if (fid > 0) H5Fclose(fid);
    delete [] xbuf;
    delete [] vbuf;
}


// Write rigid trajectories
void trajectory() 
{
    char fn[512];
    sprintf(fn, "%s%s%6.6d%s%6.6d%s", 
    		out_dir.c_str(), "rigid_traj.", 
		lt_min, "_to_", lt_max, ".dat");
    printf("Write trajectory to %s\n", fn);

    FILE *file = fopen(fn, "w");
    fprintf(file, "variables = t, x, y, z, u, v, w\n");

    MArray<double,2> x(nrec,3), v(nrec,3);

    for (int irigid = 0;  irigid < nrigid; irigid++) {
	for (int irec = 0; irec < nrec; irec++) {
	    FOR_J3 {
		x(irec,j) = xrigid(irec,irigid,j);
		v(irec,j) = vrigid(irec,irigid,j);
	    }

	    if (irec > 0) {
	        double xx = x(irec,0) - x(irec-1,0);
		xx = nearbyint(xx*ewald::iL[0])*ewald::L[0];
	        x(irec,0) = x(irec,0) - xx;
	    }
	}

	fprintf(file, "zone\n");
	for (int irec = 0; irec < nrec; irec+=10) {
	    fprintf(file, " %12.5E", Ts*step(irec));
	    fprintf(file, " %12.5E %12.5E %12.5E", 
			x(irec,0), x(irec,1), x(irec,2));
	    fprintf(file, " %12.5E %12.5E %12.5E\n", 
			v(irec,0), v(irec,1), v(irec,2));
	}
    }
    fclose(file);
}


// Density profile
void density()
{
    const double Lz = ewald::L[2];

    int nz;
    printf("Calculate density profile, nz = ? ");
//    cin >> nz;
    nz = 168;
    printf("Calculate density profile, nz = %d\n", nz);

    char fn[256];
    sprintf(fn, "%s%s%6.6d%s%6.6d%s",
		out_dir.c_str(), "rigid_dens.", 
		lt_min, "_to_", lt_max, ".dat");
    FILE *file = fopen(fn, "w");
    printf("Write rigid number desnity to %s\n", fn);
    fprintf(file, "variables = z, dens\n");

    const int LT_CHUNK = 1E5;
    MArray<double,1> dens(nz);

    dens = 0.0;
    int cnt = 0;
    for (int lt0 = step(0), irec = 0; irec < nrec; irec++) {
        int lt = step(irec);

        if ((lt/LT_CHUNK != lt0/LT_CHUNK || irec == nrec-1)  && cnt > 0) {

	    for (int iz = 0; iz < nz; iz++) {
	        dens(iz) /= cnt;
		dens(iz) /= double(nrigid)/nz;;
	    }


//	    // Force symmetry
//	    for (int iz = 0; iz < nz/2; iz++) {
//	        double tmp = 0.5*(dens(iz) + dens(nz-1-iz));
//		dens(iz) = dens(nz-1-iz) = tmp;
//	    }

	    fprintf(file, "ZONE T=\"t =%12.5E\"\n", 0.5*(lt0 + lt)*Ts);
	    for (int iz = 0; iz < nz; iz++) {
	        fprintf(file, "%12.5E %12.5E\n", (iz+0.5)*Lz/nz, dens(iz));
	    }

	    lt0 = lt;
	    dens = 0.0;
	    cnt = 0;
	}


	for (int irigid = 0; irigid < nrigid; irigid++) {
	    int iz = (int)floor(nz*xrigid(irec,irigid,2)/Lz);
	    dens(iz)++;
	}
	cnt++;
    }

    fclose(file);
} 


/* Mean velocity and velocity fluctuatio */
void velProf()
{
    const double Lz = ewald::L[2];

    int nz;
//    printf("Calculate velocity profile, nz = ? ");
//    cin >> nz;
    nz = 24;
    printf("Calculate velocity profile, nz = %d\n", nz);

    MArray<double,1> umean(nz), zmean(nz), wmean(nz), wrms(nz);
    MArray<int,1> cnt(nz);

    umean = 0.0;
    zmean = 0.0;
    wmean = 0.0;
    wrms = 0.0;
    cnt = 0;

    for (int irec = 0; irec < nrec; irec++)
    for (int irigid = 0; irigid < nrigid; irigid++) {
        double ztmp = xrigid(irec,irigid,2);
        int iz = (int)floor(nz*ztmp/Lz);
	assert(iz >= 0 && iz < nz);

	zmean(iz) += ztmp;
	umean(iz) += vrigid(irec,irigid,0);
	wmean(iz) += vrigid(irec,irigid,2);
	wrms(iz) += square( vrigid(irec,irigid,2) );
	cnt(iz)++;
    }

    for (int iz = 0; iz < nz; iz++) {
        zmean(iz) /= cnt(iz) + EPS;
	umean(iz) /= cnt(iz) + EPS;
	wmean(iz) /= cnt(iz) + EPS;
	wrms(iz) = sqrt( wrms(iz)/(cnt(iz) + EPS) );
    }

    char fn[512];
    sprintf(fn, "%s%s%6.6d%s%6.6d%s", 
    		out_dir.c_str(), "rigid_vel.",
		lt_min, "_to_", lt_max, ".dat");
    printf("Write velocity profile to %s\n", fn);

    FILE *file = fopen(fn, "w");
    fprintf(file, "VARIABLES = z, u, w, wrms\n");
    fprintf(file, "ZONE\n");
    for (int iz = 0; iz < nz; iz++) {
        if (cnt(iz) == 0) continue;
        fprintf(file, "%12.5E %12.5E %12.5E %12.5E\n", 
			zmean(iz), umean(iz), wmean(iz), wrms(iz));
    }
    fclose(file);
}


void correlation()
{
    char fn[512];
    FILE *file;

    double max_corr_time = param::getDoubleValue("CORRELATION");
    int ncorr = (int)floor(max_corr_time/(Ts*lt_step));
    ncorr = min(ncorr, nrec-10);
    printf("Max correlation time = %7.2E\n", ncorr*lt_step*Ts);

    int nz;
//    printf("Atuocorrelation, nz = ? ");
//    cin >> nz;
    nz = 12;
    printf("Atuocorrelation, nz = %d\n", nz);
    double dz = ewald::L[2]/nz;

    MArray<double,2> dz2(nz+1,ncorr), ww(nz+1,ncorr), wght(nz+1,ncorr);
    MArray<double,1> dz2_sum(ncorr), wght_sum(ncorr);

    ww = 0.0;
    dz2 = 0.0;
    wght = 0.0;
    dz2_sum = 0.0;
    wght_sum = 0.0;


    for (int irigid = 0; irigid < nrigid; irigid++) {
	// Timing
	if (irigid%5 == 0) {
	    static time_t last_t = 0;
	    time_t t = time(NULL);

	    if (irigid > 0) {
		printf("%5dth rigid processed", irigid);
		printf("    time cost = %d sec\n", (int)(t-last_t));
	    }
	    last_t = t;
        }

        MArray<double,1> z(nrec), w(nrec);
	cblas_dcopy(nrec, &xrigid(0,irigid,2), xrigid.stride(0), z.data(), 1);
	cblas_dcopy(nrec, &vrigid(0,irigid,2), vrigid.stride(0), w.data(), 1);

	for (int i0 = 0; i0 < nrec; i0 += 10) {
	    double z0 = z(i0);
	    double w0 = w(i0);
	    int iz = (int)nearbyint(z0/dz);

	    for (int j = 0; j < min(nrec,ncorr); j++) {
	        int i1 = i0 + j;
		if (i1 >= nrec) continue;

		double z1 = z(i1);
		double w1 = w(i1);

		dz2(iz,j) += square(z1 - z0);
		ww(iz,j) += w0*w1;
		wght(iz,j)++;

		if (z0 > 0.24*ewald::L[2] && z0 < 0.76*ewald::L[2]) {
		    dz2_sum(j) += square(z1 - z0);
		    wght_sum(j)++;
		}
	    }
        }
    }


    // Binned average and force symmetry
    for (int iz = 0; iz < nz; iz++) {
	for (int j = 0; j < ncorr; j++) {
	    dz2(iz,j) /= wght(iz,j) + EPS;
	    ww(iz,j) /= wght(iz,j) + EPS;
	}
    }

    for (int iz = 0; iz < nz/2; iz++) {
        int iz1 = nz - iz;

	for (int j = 0; j < ncorr; j++) {
	    double tmp;

	    tmp = ww(iz,j) + ww(iz1,j);
	    ww(iz,j) = ww(iz1,j) = 0.5*tmp;

	    tmp = dz2(iz,j) + dz2(iz1,j);
	    dz2(iz,j) = dz2(iz1,j) = 0.5*tmp;
	}
    }


    sprintf(fn, "%s%s", out_dir.c_str(), "rigid_dz2_mean.dat");
    file = fopen(fn, "w");
    printf("Write averaged displacement to %s\n", fn);
    fprintf(file, "VARIABLES = t, dz2\n");
    for (int j = 0; j < ncorr; j++) {
	fprintf(file, "%12.5E %12.5E\n", j*lt_step*Ts, dz2_sum(j)/(wght_sum(j)+1.E-10));
    }
    fclose(file);


    // Output
    sprintf(fn, "%s%s", out_dir.c_str(), "rigid_corr.dat");
    file = fopen(fn, "w");
    printf("Write correlation to %s\n", fn);
    fprintf(file, "VARIABLES = t, ww, dz2\n");

    for (int iz = 0; iz <= nz/2; iz++) {
	fprintf(file, "ZONE T=\"z = %6.3F\"\n", iz*ewald::L[2]/nz );

	for (int j = 0; j < ncorr; j++) {
	    fprintf(file, "%12.5E %12.5E %12.5E\n", 
	    		Ts*j*lt_step,
			ww(iz,j)/(ww(iz,0) + EPS),
			dz2(iz,j) );
	}
    }
    fclose(file);
}


void margination()
{
    char fn[512];
    sprintf(fn, "%s%s", out_dir.c_str(), "margination.dat");
    FILE *file = fopen(fn, "w");
    printf("Write margination to %s\n", fn);
    fprintf(file, "VARIABLES = t, z\n");

    for (int irec = 0; irec < nrec; irec++) {
        double mom = 0.0;
	for (int irigid = 0; irigid < nrigid; irigid++) {
	    double ztmp = xrigid(irec,irigid,2);
	    mom += square(ztmp - 0.5*ewald::L[2]);
	}

	mom = mom/nrigid;
	fprintf(file, "%12.5E %12.5E\n", Ts*irec*lt_step, sqrt(mom));
    }
    fclose(file);
}
