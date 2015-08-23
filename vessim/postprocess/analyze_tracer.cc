// Analyze tracers
#include "cxxheaders.h"
#include "stokesys.h"
#include "ewald.h"
#include "mathfunc.h"
#include "param.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "correlation.h"
#include "miscUtils.h"
#include <time.h>
//#include <sys/types.h>

StokeSys stks;
string in_dir("D15sean/ht10/ca05/");
string out_dir("D15sean/anatrace/");
int lt_min, lt_max, lt_step;
double Ts;

int nrec, ntrac;
MArray<double,3> xtrac, vtrac;	// (irec, itrac, xyz)
MArray<int,1> step;

void readTracer();
void velProf();
void trajectory();
void density();
void correlation();

int main()
{
    // Parameters
    param::readFile("Input/channel.in");
    param::readFile("Input/analyze_tracer.in");

    Ts = param::getDoubleValue("Ts");
    lt_min = 3000;//param::getIntValue("LT_MIN");
    lt_max = 30000;//param::getIntValue("LT_MAX");
    lt_step = 50;//param::getIntValue("LT_STEP");

    if (param::exist("IN_DIR")) in_dir = param::getStringValue("IN_DIR") + "/";
    if (param::exist("OUT_DIR")) out_dir = param::getStringValue("OUT_DIR") + "/";

    string restart_file = "D15sean/ht10/ca05/restart000000.dat";//param::getStringValue("RESTART");
    stks.readRestart(restart_file.c_str());
    // Only need to set up Ewald parameters, not the whole ewald::initialHook
    ewald::alpha = 1.0;	// just need some positive number here
    ewald::tol = 0.99;
    ewald::set_prms();

    // Alloc memory
    int nrecMax = (lt_max - lt_min)/lt_step + 1;
    ntrac = stks.ntrac;
    xtrac.resize(nrecMax,ntrac,3);
    vtrac.resize(nrecMax,ntrac,3);
    step.resize(nrecMax);

    // Read tracers
    readTracer();

    // Create directory if it doesnot exist yet
    miscUtils::mkdir(out_dir.c_str());

    if (param::exist("TRAJECTORY")) trajectory();
    if (param::exist("DENSITY")) density();
    if (param::exist("VEL_PROFILE")) velProf();
    if (param::exist("CORRELATION")) correlation();
}


void readTracer()
{
    const int LT_CHUNK = 10000;
    char fn[256], token[256];
    hid_t fid = -1;
    double *xbuf = new double[3*ntrac];
    double *vbuf = new double[3*ntrac];

    // Turn off error message
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    nrec = 0;
    for (int lt = lt_min; lt <= lt_max; lt += lt_step) {
	static int lt_file_old = -1;
	int lt_file = (lt/LT_CHUNK)*LT_CHUNK;

	// Open a new file when necessary
	if (lt_file > lt_file_old) {
	    if (fid > 0) H5Fclose(fid);

	    sprintf(fn, "%s%s%6.6d%s", in_dir.c_str(), "tracer", lt_file, ".dat");
	    printf("\t%s\n", fn);

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
	m_dcopy(3*ntrac, xbuf, &xtrac(nrec,0,0));
	m_dcopy(3*ntrac, vbuf, &vtrac(nrec,0,0));
	nrec++;
    } // lt

    printf("Number of records = %d\n", nrec);

    if (fid > 0) H5Fclose(fid);
    delete [] xbuf;
    delete [] vbuf;
}


void trajectory()
{
    double z0_min = param::getDoubleValue("TRAJECTORY", 0);
    double z0_max = param::getDoubleValue("TRAJECTORY", 1);

    char fn[512];
    sprintf(fn, "%s%s", out_dir.c_str(), "tracer_traj.dat");
    printf("Write trajectory to %s\n", fn);

    FILE *file = fopen(fn, "w");
    fprintf(file, "VARIABLES = t, x, y, z, u, v, w\n");

    MArray<double,2> x(nrec,3), v(nrec,3);

    for (int cnt = 0, itrac = 0; itrac < ntrac; itrac++) {
	for (int irec = 0; irec < nrec; irec++) {
	    FOR_J3 {
		x(irec,j) = xtrac(irec,itrac,j);
		v(irec,j) = vtrac(irec,itrac,j);
	    }

	    if (irec > 0) {
	        double xx = x(irec,0) - x(irec-1,0);
		xx = nearbyint(xx*ewald::iL[0])*ewald::L[0];
	        x(irec,0) = x(irec,0) - xx;
	    }
	}

	if (x(0,2) < z0_min || x(0,2) > z0_max) continue;

	fprintf(file, "zone\n");
	for (int irec = 0; irec < nrec; irec++) {
	    fprintf(file, " %15.5E", Ts*step(irec));
	    fprintf(file, " %15.5E %15.5E %15.5E", 
			x(irec,0), x(irec,1), x(irec,2));
	    fprintf(file, " %15.5E %15.5E %15.5E\n", 
			v(irec,0), v(irec,1), v(irec,2));
	}

	cnt++;
	if (cnt >= 100) break;
    }

    fclose(file);
}


void velProf()
{
    int nz = 24;
/*
    printf("Velocity profile, nz = %d ");
    cin >> nz;
*/
    printf("Velocity profile, nz = %d\n", nz);

    const double dz = ewald::L[2]/nz;
    const double idz = 1.0/dz;

    MArray<double,1> zmean(nz), umean(nz), wmean(nz), wrms(nz);
    MArray<int,1> cnt(nz);

    zmean = 0.0;
    umean = 0.0;
    wmean = 0.0;
    wrms = 0.0;
    cnt = 0;
    int countall = 0;
    double umeanall = 0;
    double umeantemp = 0; //Guard against roundoff error
    for (int irec = 0; irec < nrec; irec++)
    for (int i = 0; i < ntrac; i++) {
	int iz = (int)floor(xtrac(irec,i,2)*idz);

	if (iz >= 0 && iz < nz) {
	    zmean(iz) += xtrac(irec,i,2);

	    umean(iz) += vtrac(irec,i,0); 
	    umeantemp += vtrac(irec,i,0);

	    wmean(iz) += vtrac(irec,i,2);
	    wrms(iz) += square(vtrac(irec,i,2));
	    cnt(iz)++;
	    countall++;
        }
	umeanall += umeantemp;
	umeantemp = 0;
    }
    umeanall = umeanall / double(countall);
    cout << "umeanall = " << umeanall << endl;

    for (int iz = 0; iz < nz; iz++) {
        if (cnt(iz) > 0) {
	    zmean(iz) /= cnt(iz);
	    umean(iz) /= cnt(iz);
	    wmean(iz) /= cnt(iz);
	    wrms(iz) = sqrt(wrms(iz)/cnt(iz));
        }
    }

    // Write to file
    char fn[512];
    sprintf(fn, "%s%s", out_dir.c_str(), "tracer_vel.dat");
    printf("Write vel profile to %s\n", fn);

    FILE *file = fopen(fn, "w");
    fprintf(file, "variables = z, u, w, wrms\n");
    for (int iz = 0; iz < nz; iz++) {
	fprintf(file, "%12.5E %12.5E %12.5E %12.5E\n", zmean(iz), umean(iz), wmean(iz), wrms(iz));
    }

    double weightmeanall = 0;
    for (int iz = 0; iz < nz; iz++) 
      {weightmeanall += umean(iz);}
    
    weightmeanall = weightmeanall / nz;
    cout << "Height-normalized umeanall = " << weightmeanall << endl;

    fclose(file);
}


// Density profile
void density()
{
    const double Lz = ewald::L[2];

/*
    int nz;
    printf("Density profile, nz = ? ");
    cin >> nz;
*/
    int nz = 24;
    printf("Density profile, nz = %d\n", nz);

    char fn[256];
    sprintf(fn, "%s%s%6.6d%s%6.6d%s",
		out_dir.c_str(), "tracer_dens.", 
		lt_min, "_to_", lt_max, ".dat");
    printf("Write tracer number desnity to %s\n", fn);

    FILE *file = fopen(fn, "w");
    fprintf(file, "variables = z, dens\n");

    const int LT_CHUNK = 10000;
    MArray<double,1> dens(nz);

    dens = 0.0;
    int cnt = 0;

    for (int lt0 = step(0), irec = 0; irec < nrec; irec++) {
        int lt = step(irec);

        if ((lt/LT_CHUNK != lt0/LT_CHUNK || irec == nrec-1)  && cnt > 0) {
	    for (int iz = 0; iz < nz; iz++) {
	        dens(iz) /= cnt*double(ntrac)/nz;;
	    }

	    for (int iz = 0; iz <= nz/2; iz++) {
	        double tmp = 0.5*(dens(iz) + dens(nz-1-iz));
		dens(iz) = dens(nz-1-iz) = tmp;
	    }


	    fprintf(file, "zone T=\"t=%12.5E\"\n", 0.5*(lt0 + lt)*Ts);
	    for (int iz = 0; iz < nz; iz++) {
	        fprintf(file, "%12.5E %12.5E\n", (iz+0.5)*Lz/nz, dens(iz));
	    }

	    lt0 = lt;
	    dens = 0.0;
	    cnt = 0;
	}


	for (int itrac = 0; itrac < ntrac; itrac++) {
	    int iz = (int)floor(nz*xtrac(irec,itrac,2)/Lz);
	    dens(iz)++;
	}
	cnt++;
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

    // Bin the channel into nz slices
/*
    int nz;
    printf("nz = ? ");
    cin >> nz;
*/
    int nz = 12;
    printf("nz = %d\n", nz);

    double dz = ewald::L[2]/nz;
    double idz = 1.0/dz;

    MArray<double,2> ww(nz+1,ncorr), dz2(nz+1,ncorr), wght(nz+1,ncorr);
    MArray<double,1> dz2_sum(ncorr), wght_sum(ncorr);

    ww = 0.0;
    dz2 = 0.0;
    wght = 0.0;
    dz2_sum = 0.0;
    wght_sum = 0.0;

    for (int itrac = 0; itrac < ntrac; itrac++) {
	// Timing
	if (itrac%1000 == 0) {
	    static time_t last_t = 0;
	    time_t t = time(NULL);

	    if (itrac > 0) {
		printf("%5dth tracer point processed", itrac);
		printf("    time cost = %3dsec\n", (int)(t-last_t));
	    }

	    last_t = t;
        }

        MArray<double,1> z(nrec), w(nrec);
	cblas_dcopy(nrec, &xtrac(0,itrac,2), xtrac.stride(0), &z(0), 1);
	cblas_dcopy(nrec, &vtrac(0,itrac,2), vtrac.stride(0), &w(0), 1);

	for (int i0 = 0; i0 < nrec; i0 += 10) {
	    double z0 = z(i0);
	    double w0 = w(i0);
	    int iz = (int)nearbyint(z0/dz);

	    for (int j = 0; j < min(nrec, ncorr); j++) {
	        int i1 = i0 + j;
		if (i1 >= nrec) continue;
		double z1 = z(i1);
		double w1 = w(i1);

		dz2(iz,j) += square(z1 - z0);
		ww(iz,j) += w0*w1;
		wght(iz,j)++;

                if (z0 > 1./4*ewald::L[2] && z0 < 3./4*ewald::L[2]) {
		    dz2_sum(j) += square(z1 - z0);
	            wght_sum(j)++;
	        }
	    }
	} // i0
    }  // itrac


    // Mean square displacement  over the whole channel width
    sprintf(fn, "%s%s", out_dir.c_str(), "tracer_dz2_mean.dat");
    file = fopen(fn, "w");
    printf("Write spatially averaged displacement to %s\n", fn);
    fprintf(file, "VARIABLES = t, dz2\n");
    for (int j = 0; j < ncorr; j++) {
	fprintf(file, "%12.5E %12.5E\n", j*lt_step*Ts, dz2_sum(j)/wght_sum(j));
    }
    fclose(file);

    // Compute binned average
    for (int iz = 0; iz <= nz; iz++) {
	for (int j = 0; j < ncorr; j++) {
	    if (wght(iz,j) > 1.E-5) {
		dz2(iz,j) /= wght(iz,j);
		ww(iz,j) /= wght(iz,j);
	    }
        }
    }

    // Force symmetry
    for (int iz = 0; iz < nz/2; iz++) {
	int iz1 = nz - iz;

	for (int j = 0; j < ncorr; j++) {
	    double tmp;

	    tmp = 0.5*(dz2(iz,j) + dz2(iz1,j));
	    dz2(iz,j) = dz2(iz1,j) = tmp;

	    tmp = 0.5*(ww(iz,j) + ww(iz1,j));
	    ww(iz,j) = ww(iz1,j) = tmp;
        }
    }

    sprintf(fn, "%s%s", out_dir.c_str(), "tracer_corr.dat");
    file = fopen(fn, "w");
    printf("Write correlation to %s\n", fn);
    fprintf(file, "VARIABLES = t, ww, tau, dz2\n");
    for (int iz = 0; iz <= nz/2; iz++) {
	fprintf(file, "ZONE T=\"z = %6.3F\"\n", iz*dz );

	double GK = 0;

	for (int j = 0; j < ncorr; j++) {
	    double t = Ts*j*lt_step;
	    if (j > 0) {
	        GK += 0.5*lt_step*Ts*(ww(iz,j-1) + ww(iz,j))/(ww(iz,0) + 1.E-10);
            }
	    
	    fprintf(file, "%12.5E %12.5E %12.5E %12.5E\n", 
	    		t, 
			ww(iz,j)/(ww(iz,0) + 1.E-10), 
			GK,
			dz2(iz,j) );
	}
    }
    fclose(file);
} 
