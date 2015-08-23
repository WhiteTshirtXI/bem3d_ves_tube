// Read probe files calculate mean velocity profile
#include "cxxheaders.h"
#include "mathfunc.h"
#include "param.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "miscUtils.h"

double Ts;
int lt_min, lt_max, lt_step;
string in_dir("D15sean/ht10/ca1/");
string out_dir("D15sean/anaprobe/");

// Mesh of probe points
int nx, ny, nz;
vector<double> xmesh, ymesh, zmesh;

int nrec, nprb;
MArray<double,2> xprb(0,0);
MArray<double,3> vprb(0,0,0);	// (irec, iprb, xyz)
MArray<int,1> step;

void readProbe();
void velProf();

int main()
{
    param::readFile("Input/channel.in");
    param::readFile("Input/analyze_probe.in");

    Ts = param::getDoubleValue("Ts");
    lt_min = 3000;//param::getIntValue("LT_MIN");
    lt_max = 30000;//param::getIntValue("LT_MAX");
    lt_step = 50;//param::getIntValue("LT_STEP");

    if (param::exist("IN_DIR")) in_dir = param::getStringValue("IN_DIR") + "/";
    if (param::exist("OUT_DIR")) out_dir = param::getStringValue("OUT_DIR") + "/";

    readProbe();

    // Create output dir if it doesnot exist
    miscUtils::mkdir(out_dir.c_str());

    velProf();
}


/* Read the probe file */
void readProbe()
{
    const int LT_CHUNK = 10000;
    char fn[512], token[512];
    hsize_t dims[2];
    H5T_class_t class_id;
    size_t type_size;

    hid_t fid = -1;
    // Disable error message
    H5Eset_auto(H5E_DEFAULT, NULL, NULL);

    nrec = 0;
    for (int lt = lt_min; lt <= lt_max; lt += lt_step) {
	static int lt_file_old = -1;
	int lt_file = (lt/LT_CHUNK)*LT_CHUNK;

	// Open a new file when necessary
	if (lt_file > lt_file_old) {
	    if (fid > 0) H5Fclose(fid);

	    sprintf(fn, "%s%s%6.6d%s", in_dir.c_str(), "probe", lt_file, ".dat");
	    printf("\t%s\n", fn);

	    fid = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);
	    lt_file_old = lt_file;
	}
	if (fid <= 0) continue;

	// Determine the number of probes
	if (xprb.size() == 0) {

	    sprintf(token, "X");
	    if (H5LTfind_dataset(fid, token)) {
		H5LTget_dataset_info(fid, token, dims, &class_id, &type_size);
		nprb = dims[0];

		xprb.resize(nprb,3);
		H5LTread_dataset_double(fid, token, xprb.data());
	    }
	}

	sprintf(token, "V%6.6d", lt);
	if (H5LTfind_dataset(fid, token)) {
	    if (vprb.size() == 0) {
		H5LTget_dataset_info(fid, token, dims, &class_id, &type_size);
		nprb = dims[0];
		int nrec_max = (lt_max - lt_min)/lt_step + 1;
		vprb.resize(nrec_max, nprb, 3);
		step.resize(nrec_max);
	    }

	    H5LTread_dataset_double(fid, token, &vprb(nrec,0,0));
	    step(nrec) = lt;
	    nrec++;
	}
    } // lt
    if (fid > 0) H5Fclose(fid);

    // Construct the mesh
    xmesh.resize(0);
    ymesh.resize(0);
    zmesh.resize(0);

    for (int i = 0; i < nprb; i++) {
	double x = xprb(i,0);
	double y = xprb(i,1);
	double z = xprb(i,2);
	double eps = 1.E-5;

	vector<double>::iterator it;

	it = std::upper_bound(xmesh.begin(), xmesh.end(), x-eps);
	if (it == xmesh.end() || *it > x+eps) xmesh.insert(it, x);

	it = std::upper_bound(ymesh.begin(), ymesh.end(), y-eps);
	if (it == ymesh.end() || *it > y+eps) ymesh.insert(it, y);

	it = std::upper_bound(zmesh.begin(), zmesh.end(), z-eps);
	if (it == zmesh.end() || *it > z+eps) zmesh.insert(it, z);
    } // i

    nx = xmesh.size();
    ny = ymesh.size();
    nz = zmesh.size();
    printf("Number of probes = %d\n", nprb);
    printf("Number of records = %d\n", nrec);

    printf("Probe mesh = %d x %d x %d\n", nx, ny, nz);
}


void velProf()
{
    MArray<double,1> umean(nz), wrms(nz);
    MArray<int,1> cnt(nz);

    umean = 0.0;
    wrms = 0.0;
    cnt = 0;

    for (int iprb = 0; iprb < nprb; iprb++) {
	int iz;
	for (iz = 0; iz < nz; iz++) {
	    if ( fabs(zmesh[iz] - xprb(iprb,2)) < 1.E-10 ) break;
	}

	for (int irec = 0; irec < nrec; irec++) {
	    umean(iz) += vprb(irec,iprb,0);
	    wrms(iz) += square(vprb(irec,iprb,2));
	}
	cnt(iz) += nrec;
    }

    for (int iz = 0; iz < nz; iz++) {
        umean(iz) /= cnt(iz);
	wrms(iz) = sqrt(wrms(iz)/cnt(iz));
    }

    char fn[512];
    sprintf(fn, "%s%s", out_dir.c_str(), "probe_vel.dat");
    printf("Write velocity profile to %s\n", fn);

    FILE *file = fopen(fn, "w");
    fprintf(file, "VARIABLES = z, u, wrms\n");
    fprintf(file, "ZONE I=%d F=POINT\n", nz);
    for (int iz = 0; iz < nz; iz++) {
	fprintf(file, "%12.5E %12.5E %12.5E\n", zmesh[iz], umean(iz), wrms(iz));
    }

    double allumean = 0;
    for (int iz = 0; iz < nz; iz++) {
      allumean += umean(iz);
    }
    allumean = allumean/nz;
    cout << "Overall mean velocity = " << allumean << endl;

    fclose(file);
}
