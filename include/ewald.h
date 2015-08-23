#ifndef EWALD_H
#define EWALD_H

#include "classdec.h"
#include "rfftw_mpi.h"
#include "petscmat.h"

namespace ewald
{
    // Bounary condition type
    //   REC -- rectangular domain
    //   SHEAR -- simple shear
    //   STRAIN -- planar strain
    // Note: 
    //   For shear flow,  (u, v, w) = gamma * (z, 0, 0)
    //   For strain flow, (u, v, w) = gamma * (x, 0, -z)
    enum BC_TYPE { RECT = 0, SHEAR = 1, STRAIN = 2 };
}

namespace ewald
{
    extern BC_TYPE bctype;		// boundary condition type
    extern double axis[3][3]; 		// unit cell axis
    extern double rcp_axis[3][3];	// reciprocal axis
    extern double strain;		// the total strain

    extern double L[3];			// unit cell size when strain = 0
    extern double iL[3];		// iL = 1/L
    extern double vol, iVol;		// volume, iVol = 1/vol

    extern double alpha;		// Ewald sum alpha
    extern double tol;			// Ewald sum error tol
    
    void init_comm();
    void finalize_comm();
    void set_prms();

    void updateThreeAxes();
    void phys_to_lattice_coord(const double *, double *);
    void lattice_to_phys_coord(const double *, double *);

    void calcSum(const DoubArray2d &x, DoubArray2d &v,
    		const DoubArray2d &y, 
		const DoubArray2d &f, const DoubArray2d &g, const DoubArray2d &nrml,
		double tol);

    void calcPhysSum(const DoubArray2d &x, DoubArray2d &v,
    		const DoubArray2d &y, 
		const DoubArray2d &f, const DoubArray2d &g, const DoubArray2d &nrml,
		double tol);

    void calcFourSum(const DoubArray2d &x, DoubArray2d &v,
    		const DoubArray2d &y, 
		const DoubArray2d &f, const DoubArray2d &g, const DoubArray2d &nrml,
		double tol);


    void initialHook();
    void finalHook();

    void to_CloseBy(const double *xref, double *x);

    // Timing
    extern bool timing;
    void getWallTime(double &, double &);
}


namespace ewald {
namespace phys {
    extern bool active;
    extern MPI_Comm comm;
    extern int comm_size, comm_rank;

    extern double rc;

    void coeff_exact(double, double *, double *, double *);
    void coeff(double, double *, double *, double *);

    void int_on_tri_regular(const double *, const Tri &,
            double (*lhs1)[3][3]=NULL, double(*lhs2)[3][3]=NULL,
            double c1=0.0, double c2=0.0, double* rhs=NULL,
            const char *qrule="TRI_3");

    void int_on_tri_singular(const double *, const Tri &,
            double (*lhs1)[3][3]=NULL, double(*lhs2)[3][3]=NULL,
            double c1=0.0, double c2=0.0, double* rhs=NULL,
            int nqp=6);

    void int_on_tri_nearsingular(const double *, const Tri &,
            double (*lhs1)[3][3]=NULL, double(*lhs2)[3][3]=NULL,
            double c1=0.0, double c2=0.0, double* rhs=NULL,
            int nqp=6);

    void addSurfInt(NbrList&, double, double, double*);
    void calcSurfIntMat(NbrList&, Mat, Mat);

    // Timing
    void resetTiming();
    void startTiming();
    void addWallTime();
} }


namespace ewald {
namespace four
{
    extern bool active;
    extern MPI_Comm comm;
    extern int comm_size, comm_rank;

    extern double xmin, xmax, xminBuf, xmaxBuf;

    extern const int PB;	// (B-spline order) + 1
    extern int Nb[3];

    extern int iBgn[3], iEnd[3];
    extern int qBgn[3], qEnd[3];

    extern bool _f_flag, _g_flag;
    extern MArray<double,4>  _f, _g, _v;
    extern MArray<double,1> bb0, bb1, bb2;

    extern rfftwnd_mpi_plan _fplan, _bplan;

    bool point_active(const double *);
    bool tri_active(const double (*)[3]);

    void setSources(const vector<Tri*> &, vector<Tri*> &);
    void setTargets(const vector<Point*> &, vector<Point*> &);
    void getBlkNum(const double *, double *);

    void initPme();
    void clear_source();
    void add_source(int, const double *, double, const double *,
            double, const double *, const double *);
    void add_source(double, double, const vector<Tri*> &);

    void add_interp_vel(int, const double *, double *);
    void add_interp_vel(const vector<Point*> &, double *);
    void transform();

    // Timing
    void resetTiming();
    void startTiming();
    void addWallTime();
} }


// Source list
class SrcList
{
public:
    vector<Tri*> faces;

    double lb[3], ub[3];  // domain boundaries
    bool cyclic[3];       // periodicity flag

    // Verlet cell list
    // The box contains all the source triangles and a buffer zone
    // of rc width.  A target point must have Verlet cell index in
    // [0,nblk) to have any interaction with sources.
    int nblk[3];          // number of blocks
    double blkSz[3], iblkSz[3];   // block size
    MArray<int,3> hoc;
    MArray<int,1> next;

    SrcList() { }
    ~SrcList() { }

    int numFaces() { return faces.size(); };

    void build(const vector<Tri*> &allFaces);

    void buildBlks();
    void buildBlks_rect();
    void buildBlks_general();

    void getBlkNum(const double *, int &, int &, int &);
    void getBlkNum_rect(const double *, int &, int &, int &);
    void getBlkNum_general(const double *, int &, int &, int &);
};


// Neighbor list for short range interaction
class NbrList
{
public:
    vector<Point*> verts;
    vector<Tri*> faces;

    vector<double> dists;
    vector<int> firstNbr;

    NbrList() { };
    ~NbrList() { };

    void build(const vector<Point*> &, SrcList &,
                    bool NO_SELF=false, bool ONLY_SELF=false);
    void build_rect(const vector<Point*> &, SrcList &,
                    bool NO_SELF=false, bool ONLY_SELF=false);
    void build_general(const vector<Point*> &, SrcList &,
                    bool NO_SELF=false, bool ONLY_SELF=false);

    void build(const vector<Point*> &, const vector<Tri*> &);
    void findNumNbrPoints(int*);
};

#endif
