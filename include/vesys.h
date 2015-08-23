#ifndef VESYS_H
#define VESYS_H

#include "classdec.h"
#include "vesicle.h"
#include "ewald.h"
#include "mypetsc.h"

class VeSys {
public:
    vector<Vesicle> vesicles;

    vector<Mesh*> meshList;
    vector<Point*> vertList;
    vector<Tri*> faceList;

    // For Ewald physical sum
    SrcList slist_phys;
    NbrList nlist_phys;

    // For Ewald Fourier sum
    vector<Tri*> slist_four;
    vector<Point*> tlist_four;

    Mat matSL, matDL;

    int Nt0, Nt, lt;
    double Ts, time;

    double shRate;		// shear rate
    double viscRat;		// viscosity ratio

    VeSys();
    ~VeSys();

    // Background shear flow velocity
    virtual void calcBkgVel(const double *x, double *v) const {
        v[0] = shRate*(x[2] - 0.5*ewald::L[2]);
	v[1] = 0.0;
	v[2] = 0.0;
    }

    int numVesicles() const { return vesicles.size(); }

    void initialHook();

    void timeInt();

    void cellNoContact(double);

    void domainDecomp();
    void updateSourceGeometry();
    void updateNeighborList();

    void calcBieMat();

    void rebox();
    void syncCoord();

    void calcStress(double (*)[3]);

    void writeAll();
    void writeVesicles(const char *);
    void writeVesicleCenter(const char *);
    void writeStress(const char *, const double (*)[3]);

    void writeRestart(const char *);
    void readRestart(const char *);

    Vec rhs_vel, sol_vel;
    VecScatter scatter_vel;

    void addBoundaryIntegral(double cf, double cg, Vec);
    
    Vec rhs_lbd, sol_lbd;
    VecScatter scatter_lbd;

    void solveBie(double);
    void initBieSolver();
    void finalizeBieSolver();

    void projectVel(const DoubArray2d &, const DoubArray1d &,
    		DoubArray2d &, DoubArray1d &, bool initial_guess = false);
    void initProjectSolver();
    void finalizeProjectSolver();

    void solveVel();

    double volFrac();
    void calcMaxVelPerturb(double *dvmax);
};

#endif
