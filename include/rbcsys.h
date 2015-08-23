// Homogeneous shear flow system
#ifndef RBCSYS_H_
#define RBCSYS_H_

#include "classdec.h"
#include "redcell.h"
#include "ewald.h"
#include "mypetsc.h"
#include "hdf5.h"
#include "hdf5_hl.h"

class RbcSys
{
public:
    vector<RedCell> cells;

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

    RbcSys();
    ~RbcSys();

    // Background shear flow velocity
    virtual void calcBkgVel(const double *x, double *v) const {
        v[0] = shRate*(x[2] - 0.5*ewald::L[2]);
	v[1] = 0.0;
	v[2] = 0.0;
    } 

    int numCells() const { return cells.size(); }

    void initialHook();

    void timeInt();

    void cellNoContact(double);

    void updateGeometry();
    void updateSourceList();
    void updateTargetList();
    void updateNeighborList();
    void calcBieMat();

    void rebox();
    void syncCoord();

    void calcStress(double (*)[3]);

    void writeAll();
    void writeCells(const char *);

    void writeRestart(const char *);
    void readRestart(const char *);

    void writeCellCenter(const char *);
    void writeStress(const char *, const double (*)[3]);

    Vec rhs_vel, sol_vel;
    VecScatter scatter_vel;

    void addBoundaryIntegral(double, DoubArray2d &, double, DoubArray2d &, Vec);
    void calcResidual_Bie(Vec);

    void initBieSolver();
    void finalizeBieSolver();
    void solveBie(double);
    void verifySolution();

    double cellVolFraction();
    void calcMaxVelPerturb(double *dvmax);
    void calcAreaChange(double &dAmi, double &dAmax, double &dAmean);
};

#endif
