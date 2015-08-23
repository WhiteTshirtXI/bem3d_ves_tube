#ifndef STOKESYS_H_
#define STOKESYS_H_

#include "classdec.h"
#include "vesicle.h"
#include "rigid.h"
#include "wall.h"
#include "ewald.h"
#include "mypetsc.h"
#include "hdf5.h"
#include "hdf5_hl.h"

const int CELL = 0;
const int RIGID = 1;
const int WALL = 2;

class StokeSys
{
public:
    vector<Vesicle> cells;  //This line is now changed to type vesicles instead of cell.
    vector<Rigid> rigids;
    vector<Wall> walls;

    vector<Mesh*> meshList[3];
    vector<Point*> vertList[3];
    vector<Tri*> faceList[3];

    // For Ewald physical sum
    SrcList slist_phys[3];
    NbrList nlist_phys[3][3];

    SrcList slist_RigidSelf;
    NbrList nlist_RigidSelf; 

    // For Ewald Fourier sum
    vector<Tri*> slist_four[3];
    vector<Point*> tlist_four[3];

    int Nt, lt;
    double Ts;
    double time;

    double pgrad[3];		// pressure gradient
    double vbkg[3];		// background velcoity
    double viscRat;		// viscosity ratio

    // Probes (Eulerian) and tracers (Lagrangian)
    int nprb, ntrac;
    double (*xprb)[3], (*vprb)[3];
    double (*xtrac)[3], (*vtrac)[3];

	string outfiledirectory;

    StokeSys();
    ~StokeSys();

    int numCells() const { return cells.size(); }
    int numRigids() const { return rigids.size(); }
    int numWalls() const { return walls.size(); }

    void initialHook();

    void timeInt();

    void cellNoContact(double);
    void addCollisionForce(int, int, double);

    bool pointInsideSomeCell(const double *, int *whichCell=NULL);

    void calcFixedSTNList();
    void domainDecomp();
    void setMeshOwnership();
    void updateSourceGeometry();
    void buildNeighborList();

    void rebox();
    void syncCoord();

    void calcPressGrad(double *);
    void calcCellVolFrac(double, double, int, double *);
    void calcStress(double, double, int, 
    		double *, double *, double *, double *,
		double *, double *, double *);
    void calcStressLocal(double, double, int, 
    		double *, double *, double *, double *, 
		double *);
		
    void calcParticleDensity(double, double, int, double *);
    void calcProbeVel(int n, const double (*)[3], double (*)[3]);

    void writeAll();
    void writeCells(const char *);
    void writeRigids(const char *);
    void writeWalls(const char *);

    void writeRestart(const char *);
    void readRestart(const char *);
    void readRestart_quiet(const char *);


    void writeCellCenter(const char *);
    void writeRigidCenter(const char *);
    void writeProbe(const char *);
    void writeTracer(const char *);

    void initFlowSolver();
    void finalizeFlowSolver();
    void solveFlow();
    void verifySolution();
    
    void solveBie(double);
	void projectSol(const DoubArray2d &divTar, DoubArray1d &lbd,
		DoubArray2d &cell_v, DoubArray2d &rigid_v, DoubArray2d &wall_f,
		bool initial_guess);
		
	void initProjectSolver();
    void finalizeProjectSolver();

    double cellVolFraction();
    //void calcAreaChange(double &dAmin, double &dAmax, double &dAmean); //Cell-only function
};

#endif
