#include "cxxheaders.h"
#include "ewald.h"
#include "mathfunc.h"
#include "geom_oper.h"

/* void ewald::initialHook
 * void ewald::finalHook
 * void ewald::set_prms
 * void ewald::init_comm
 * void ewald::finalize_comm
 * void ewald::to_CloseBy
 * void ewald::getWallTime */

namespace ewald
{
    BC_TYPE bctype = RECT;
    double strain = 0.0;
    double axis[3][3], rcp_axis[3][3];

    double L[3], iL[3];
    double vol, iVol;

    // Non-sense default values to force explicit inialization
    double alpha = -1.0;
    double tol = -1.0;

    bool timing = false;

    /* Update the axes and reciprocal axes of the unit cell
     * Note:
     *   -- Needs to preset: bctype, L, strain (if bctype is SHEAR or STRAIN ) */
    void updateThreeAxes() {
	// Axis
	m_dclear(9, *axis);
	m_dclear(9, *rcp_axis);

	const double HenckyStrain = log( 0.5*(3.0 + sqrt(5.0)) );
	double s;

        switch (bctype) {
	    case RECT:
	        axis[0][0] = L[0];
	        axis[1][1] = L[1];
	        axis[2][2] = L[2];
	        break;

	    case SHEAR:
		axis[0][0] = L[0];
		axis[1][1] = L[1];
		axis[2][0] = strain*L[2] - nearbyint(strain*L[2]/L[0])*L[0];
		axis[2][2] = L[2];
	        break;

	    case STRAIN:
		// For now, can only handle the case where the undeformed shape
		// in the x-z plane is a square one

		// Init domain at zero strain
	        assert ( fabs(L[0] - L[2]) < 1.E-10*(L[0] + L[2]) );
		axis[0][0] = 1.0;
		axis[0][2] = -(sqrt(5.0) - 1)/2;
		s = 1.0/m_dnrm2(3, axis[0]);
		m_dscal(3, s*L[0], axis[0]);

		axis[1][1] = L[1];

		axis[2][0] = -axis[0][2];
		axis[2][2] = axis[0][0];

		// Hencky strain = the point where the deformed lattice
		//                 becomes identical to the original one
		s = exp( strain - floor(strain/HenckyStrain)*HenckyStrain );
		axis[0][0] *= s;
		axis[2][0] *= s;

		axis[0][2] /= s;
		axis[2][2] /= s;
	        break;

	    default:
	        printf("Invalid ewald::bctype\n");
		exit(1);
	}

	// Reciprocal axis
	FOR_I3 {
	    cross_product( axis[(i+1)%3], axis[(i+2)%3], rcp_axis[i] );
	}
	m_dscal(9, iVol, *rcp_axis);
    }


    /* Convert from physical coordinates to lattice coordinates
     * Arguments:
     *   x -- physical coordinates
     *   s -- lattice coordinates */
    void phys_to_lattice_coord(const double *x, double *s) {
        if (bctype == RECT) {
	    FOR_D3 s[d] = x[d]*iL[d];
	} 
	else {
	    s[0] = rcp_axis[0][0]*x[0] + rcp_axis[0][2]*x[2];
	    s[1] = rcp_axis[1][1]*x[1];
	    s[2] = rcp_axis[2][0]*x[0] + rcp_axis[2][2]*x[2];
	}
    }

    /* Convert from lattice coordinates to physical coordinates
     * Arguments:
     *   s -- lattice coordinates 
     *   x -- physical coordinates */
    void lattice_to_phys_coord(const double *s, double *x) {
        if (bctype == RECT) {
	    FOR_D3 x[d] = s[d]*L[d];
	} 
	else {
	    x[0] = s[0]*axis[0][0] + s[2]*axis[2][0];
	    x[1] = s[1]*axis[1][1];
	    x[2] = s[0]*axis[0][2] + s[2]*axis[2][2];
	}
    }
}


namespace ewald {
namespace phys
{
    bool active = true;
    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_size = 1; 
    int comm_rank = 0;

    double rc;

    // Timing
    bool _is_timing = false;
    double _wtime, _wtimeBgn;

    void resetTiming() 
    { 
        _wtime = 0.0; 
	_is_timing = false; 
    }

    void startTiming() 
    { 
	_wtimeBgn = MPI_Wtime(); 
	_is_timing = true; 
    }

    void addWallTime() 
    { 
        assert(_is_timing); 
	_wtime += MPI_Wtime() - _wtimeBgn; 
	_is_timing = false;
    }
}}

namespace ewald {
namespace four
{
    bool active = true;
    MPI_Comm comm = MPI_COMM_WORLD;
    int comm_size = 1; 
    int comm_rank = 0;

    double xmin, xmax, xminBuf, xmaxBuf;

    const int PB = 4;
    int Nb[3];

    int iBgn[3], iEnd[3];
    int qBgn[3], qEnd[3];

    bool _f_flag, _g_flag;
    MArray<double,4>  _f, _g, _v;
    MArray<double,1> bb0, bb1, bb2;

    rfftwnd_mpi_plan _fplan, _bplan;

    // Timing
    bool _is_timing = false;
    double _wtime, _wtimeBgn;

    void resetTiming() 
    { 
        _wtime = 0.0; 
	_is_timing = false; 
    }

    void startTiming() 
    { 
	_wtimeBgn = MPI_Wtime(); 
	_is_timing = true; 
    }

    void addWallTime() 
    { 
        assert(_is_timing); 
	_wtime += MPI_Wtime() - _wtimeBgn; 
	_is_timing = false;
    }
}};


// Ewald initial hook
void ewald::initialHook()
{
    init_comm();

    set_prms();
    updateThreeAxes();

    if (four::active) four::initPme();
}


// Ewald final hook
void ewald::finalHook()
{
    finalize_comm();
}


/* Set up Ewald parameters
 * Note:
 *   Input: L, alpha and tol */
void ewald::set_prms()
{
    int mpi_rank = 0;
    
    int mpi_inited;
    MPI_Initialized(&mpi_inited);
    if (mpi_inited) {
        MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    }

    if (mpi_rank == 0) {
	printf("Ewald sum:\n");
	printf("  L = %.3f %.3f %.3f\n", L[0], L[1], L[2]);
    }

    // Basic geometrical parameters
    FOR_I3 iL[i] = 1.0/L[i];
    vol = L[0]*L[1]*L[2];
    iVol = 1.0/vol;

    if (alpha <= 0 || tol <= 0 || tol >= 1.0) {
        if (mpi_rank == 0) {
            printf("ERROR: invalid Ewald input parameters\n");
        }
	return;
    }

    // Physical sum
    if (phys::active) {
        using phys::comm;
        using phys::comm_size;
        using phys::comm_rank;
        using phys::rc;

	double s = 1.0;
	for (int iter = 0; iter < 10; iter++) {
	    s = 0.75*sqrt(M_PI)*tol/(s*s*s + 1.5*s + 0.75/s);
	    s = sqrt(-log(s));
	}
	rc = sqrt(alpha/M_PI)*s;

	if (comm_rank == 0) {
	    printf("  rc = %.3f\n", phys::rc);
        }
    }

    // Fourier sum
    if (four::active) {
	using four::comm;
	using four::comm_size;
	using four::comm_rank;
        using four::Nb;
        using four::PB;

	double r = sqrt(-M_PI*alpha/log(tol));

	FOR_D3 {
	    Nb[d] = 2*(int)ceil(L[d]/r);
	    Nb[d] = max(Nb[d], PB);
	    if (Nb[d]%comm_size != 0) Nb[d] += comm_size - Nb[d]%comm_size;
	}

	if (comm_rank == 0) {
	    printf("  SPME mesh = %d %d %d\n", Nb[0], Nb[1], Nb[2]);
	    printf("  SPME PB = %d\n", PB);
	}
    }
}


/* Init Ewald sum MPI */
void ewald::init_comm()
{
    int mpi_inited;
    MPI_Initialized(&mpi_inited);
    if (! mpi_inited) MPI_Init(NULL, NULL);

    int mpi_size, mpi_rank;
    char proc_name[MPI_MAX_PROCESSOR_NAME];
    int len_proc_name;
    MPI_Group group_world, group_phys, group_four;
    MPI_Status stat;

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Get_processor_name(proc_name, &len_proc_name);

    if (mpi_rank == 0) {
        printf("Number of processors = %3d\n", mpi_size);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // For now, no split of communicator
    phys::active = true;
    four::active = true;

    // find the members in phys::comm and four::comm
    int flagP, flagF;
    int cntP, *ranksP, cntF, *ranksF;

    ranksP = new int[mpi_size];
    ranksF = new int[mpi_size];

    flagP = phys::active ? 1 : 0;
    flagF = four::active ? 1 : 0;

    if (mpi_rank == 0) {
        cntP = 0;
        cntF = 0;

        if (flagP > 0) ranksP[cntP++] = mpi_rank;
        if (flagF > 0) ranksF[cntF++] = mpi_rank;

        for (int i = 1; i < mpi_size; ++i) {
            MPI_Recv(&flagP, 1, MPI_INT, i, 111, MPI_COMM_WORLD, &stat);
            MPI_Recv(&flagF, 1, MPI_INT, i, 222, MPI_COMM_WORLD, &stat);

            if (flagP > 0) ranksP[cntP++] = i;
            if (flagF > 0) ranksF[cntF++] = i;
        }
    } else {
        MPI_Send(&flagP, 1, MPI_INT, 0, 111, MPI_COMM_WORLD);
        MPI_Send(&flagF, 1, MPI_INT, 0, 222, MPI_COMM_WORLD);
    }

    MPI_Bcast(&cntP, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ranksP, cntP, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(&cntF, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(ranksF, cntF, MPI_INT, 0, MPI_COMM_WORLD);

    // Create phys::comm and four::comm
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
    MPI_Group_incl(group_world, cntP, ranksP, &group_phys);
    MPI_Comm_create(MPI_COMM_WORLD, group_phys, &phys::comm);
    MPI_Group_free(&group_phys);

    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
    MPI_Group_incl(group_world, cntF, ranksF, &group_four);
    MPI_Comm_create(MPI_COMM_WORLD, group_four, &four::comm);
    MPI_Group_free(&group_four);

    // Delete temp arrays
    delete [] ranksP;
    delete [] ranksF;

    // Print info
    if (phys::active) {
        using phys::comm;
        using phys::comm_size;
        using phys::comm_rank;

        MPI_Comm_size(comm, &comm_size);
        MPI_Comm_rank(comm, &comm_rank);

        if (comm_rank == 0) {
            printf("Ewald physical sum runs on %d processors\n", comm_size);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);


    if (four::active) {
        using four::comm;
	using four::comm_size;
	using four::comm_rank;

        MPI_Comm_size(comm, &comm_size);
        MPI_Comm_rank(comm, &comm_rank);

        if (comm_rank == 0) {
            printf("Ewald Fourier  sum runs on %d processors\n", comm_size);
        } 
    }
    MPI_Barrier(MPI_COMM_WORLD);
}


// Finalize MPI communicator for Ewald sum
void ewald::finalize_comm()
{
    if (phys::active) MPI_Comm_free(&phys::comm);
    if (four::active) MPI_Comm_free(&four::comm);
}


// Move x close to xref
void ewald::to_CloseBy(const double *xref, double *x)
{
    if (bctype == RECT) {
        FOR_I3 {
	    double rr = (x[i] - xref[i])*iL[i];
	    rr = nearbyint(rr);
	    x[i] -= rr*L[i];
	}
    } 
    else {
        double xx[3], rr[3];
	FOR_I3 xx[i] = x[i] - xref[i];
	phys_to_lattice_coord(xx, rr);

	FOR_I3 rr[i] = nearbyint(rr[i]);

	lattice_to_phys_coord(rr, xx);
	FOR_I3 x[i] -= xx[i];
    }
}


/* Get wall time spent on physical and Fourier sums
 * Arguments:
 *  wtime_phys -- wall time on physical sum
 *  wtime_four -- wall time on Fourier sum */
void ewald::getWallTime(double &wtime_phys, double &wtime_four)
{
    double wtime_tmp = 0.0;
    if (phys::active && phys::comm_rank == 0) wtime_tmp = phys::_wtime;
    MPI_Allreduce(&wtime_tmp, &wtime_phys, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    wtime_tmp = 0.0;
    if (four::active && four::comm_rank == 0)  wtime_tmp = four::_wtime;
    MPI_Allreduce(&wtime_tmp, &wtime_four, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}
