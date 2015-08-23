// Homogeneous shear flow system
#ifndef STOKESYS_H_
#define STOKESYS_H_

#include "classdec.h"
#include "cell.h"
#include "rigid.h"
#include "ewald.h"
#include "mypetsc.h"
#include "hdf5.h"
#include "hdf5_hl.h"

const int CELL = 0;
const int RIGID = 1;

class ShearSys
{
public:
    vector<Cell> cells;
    vector<Rigid> rigids;

    vector<Mesh*> meshList[2];
    vector<Point*> vertList[2];
    vector<Tri*> faceList[2];

    // For Ewald physical sum
    SrcList slist_phys[2];
    NbrList nlist_phys[2][2];

    SrcList slist_RigidSelf;
    NbrList nlist_RigidSelf;

    // For Ewald Fourier sum
    vector<Tri*> slist_four[2];
    vector<Point*> tlist_four[2];

    int Nt0;
    int Nt, lt;
    double Ts;
    double time;

    double shRate;		// shear rate
    double viscRat;		// viscosity ratio

    // Probes (Eulerian) and tracers (Lagrangian)
    int nprb, ntrac;
    double (*xprb)[3], (*vprb)[3];
    double (*xtrac)[3], (*vtrac)[3];

    ShearSys();
    ~ShearSys();

    // Background shear flow velocity
    virtual void calcBkgVel(const double *x, double *v) const {
        v[0] = shRate*(x[2] - 0.5*ewald::L[2]);
	v[1] = 0.0;
	v[2] = 0.0;
    } 

    int numCells() const { return cells.size(); }
    int numRigids() const { return rigids.size(); }

    void initialHook();

    void timeInt();

    void cellNoContact(double);
    void addCollisionForce(int, int, double);
    bool pointInsideSomeCell(const double *, int *whichCell=NULL);

    void domainDecomp();
    void buildNeighborList();
    void updateSourceGeometry();

    void rebox();
    void syncCoord();

    void calcStress(double (*)[3]);
    void calcProbeVel(int n, const double (*)[3], double (*)[3]);

    void writeAll();
    void writeCells(const char *);
    void writeRigids(const char *);

    void writeRestart(const char *);
    void readRestart(const char *);

    void writeCellCenter(const char *);
    void writeRigidCenter(const char *);
    void writeProbe(const char *);
    void writeTracer(const char *);
    void writeStress(const char *, const double (*)[3]);

    void initFlowSolver();
    void finalizeFlowSolver();
    void solveFlow();

    double cellVolFraction();
    void calcMaxPerturbVel(double *dvmax, double *UTmax, double *OMGmax);
    void calcAreaChange(double &dAmi, double &dAmax, double &dAmean);
};

#endif
