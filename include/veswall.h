#ifndef VESICLE_WALL_H
#define VESICLE_WALL_H

#include "classdec.h"
#include "vesicle.h"
#include "wall.h"
#include "ewald.h"
#include "mypetsc.h"

const int VESICLE = 0;
const int WALL = 1;

class VesWall {
public:
    vector<Vesicle> vesicles;
    vector<Wall> walls;

    vector<Mesh*> meshList[2];
    vector<Point*> vertList[2];
    vector<Tri*> faceList[2];

    // For Ewald physical sum
    SrcList slist_phys[2];
    NbrList nlist_phys[2][2];

    // For Ewald Fourier sum
    vector<Tri*> slist_four[2];
    vector<Point*> tlist_four[2];

    // Interaction matries
    Mat matSL[2][2], matDL[2][2];

    // Probes (Eulerian) and tracers (Lagrangian)
    int nprb, ntrac;
    MArray<double,2> xprb, vprb;
    MArray<double,2> xtrac, vtrac;
    double (*xprb_ptr)[3];// = &xprb; 
    double (*vprb_ptr)[3];// = &vprb;
    double (*xtrac_ptr)[3];// = &xtrac;
    double (*vtrac_ptr)[3];// = &vtrac;

    string outfiledirectory;

    int Nt0, Nt, lt;
    double Ts, time;

    double vbkg[3];		// mean velocity
    double shRate;		// shear rate
    double viscRat;		// viscosity ratio

    VesWall();
    ~VesWall();

    //int numCells() const { return vesicles.size(); } //numCells is an alias for numVesicles right now for compatibility reasons
    int numVesicles() const { return vesicles.size(); }
    int numWalls() const { return walls.size(); }

    void initialHook();
    void domainDecomp();
    void updateSourceGeometry();
    void updateNeighborList();

    void calcBieMat();
    void addBoundaryIntegrals(int iblk, Vec vec, const double *, const double *);
    void addBoundaryIntegrals(const MArray<double,2> &x, Vec vec, const double *, const double *);

    void calcRhs_vesicleBie(Vec);
    void calcRhs_wallBie(Vec);

    // Vesicle BIE solver
    void solveWallBie();
    void initWallBieSolver();
    Vec getWallBieRhs();
    void solveVesicleBie(double myDt);
    void initVesicleBieSolver();
    Vec getVesicleBieRhs();
    Vec getVesicleBieSol();
    VecScatter getVesicleBieScatter();

    // Vesicle surface tension solver
    void projectVel(const DoubArray2d &vold, const DoubArray1d &divTar,
                    DoubArray2d &vnew, DoubArray1d &lbd,
		    bool initial_guess=false);
    void initProjectSolver();

    // Wall solver
    Vec getWallBieSol();
    VecScatter getWallBieScatter();

    void calcProbeVel(MArray<double,2> &x, MArray<double,2> &v);

    void timeInt();
    void rebox();
    void syncCoord();

    void writeAll();
    void writeVesicles(const char *);
    void writeWalls(const char*);
    void writeRestart(const char *);
    void writeProbe(const char*);
    void writeVesicleCenter(const char *);
    void writeStress(const char *);
    void readRestart(const char *);

    void calcStress(double (*)[3]);
    double volFraction();
    
    bool pointInsideSomeVesicle(const double *, int *whichVesicle=NULL);
};

#endif
