#ifndef MESH_H
#define MESH_H

#include "classdec.h"
#include "point.h"
#include "tri.h"
#include "edge.h"

class Mesh
{
public:
    int isActive;	// active flag
    int isPrivate;	// privately owned flag

    int Gindx;		// global index

    vector<Point> verts;
    vector<Tri> faces;
    vector<Edge> edges;

    // Area, volume, and centroid (weighted by element area)
    double area, vol, center[3];

    // pmom -- eigenvalues of moment of inertia tensor
    // paxis[i][:] -- the i-th principal axis of mom
    double pmom[3], paxis[3][3]; 

    MArray<double,1> vertArea;	// surface area assigned to vertex
    MArray<double,2> vertNrml;	// surface normal at vertex
    MArray<double,1> vertH;	// mean curvature at vertex

    Mesh();
    Mesh(const Mesh&);
    Mesh& operator=(const Mesh&);

    ~Mesh();

    int numVerts() const { return verts.size(); }
    int numFaces() const { return faces.size(); }
    int numEdges() const { return edges.size(); }

    void getCoords(double *) const;
    void getCoords(double (*)[3]) const;
    void getCoords(MArray<double,2> &) const;

    void setCoords(const double *);
    void setCoords(const double (*)[3]);
    void setCoords(const MArray<double,2> &);

    void getConnectivities(int *) const;
    void getConnectivities(int (*)[3]) const;
    void getConnectivities(MArray<int,2> &) const;

    void setConnectivities(const int *);
    void setConnectivities(const int (*)[3]);
    void setConnectivities(const MArray<int,2> &);

    void setInternalPointers();
    void connectPeriodicBoundaries(const double *L=NULL);
    void noPeriodicBoundaries();

    void updateAVC();
    void updateGeometry();

    void getCoordRange(double *, double *);
    void getCoordRange(int, double &, double &);

    void calcAreaAndVolume(double &A, double &V);

    void calcVertArea(MArray<double,1> &);
    void calcVertNormal(MArray<double,2> &, MArray<double,1> &);

    void calcMomentInertia(double *, double (*)[3]);

    void shapeFactor(double &psi, double &D);

    void calcCentroidVel(const double *, double *);
    void calcCentroidVel(const double (*)[3], double *);
    void calcCentroidVel(const MArray<double,2> &, double *);

    void calcTransRotatVel(const double *, double *, double *);
    void calcTransRotatVel(const double (*)[3], double *, double *);
    void calcTransRotatVel(const MArray<double,2> &, double *, double *);

    void calcMinMaxEdgeLen(double *, double *);
    double minEdgeLen();
    double maxEdgeLen();
    double minAngle();
    bool isInteriorPoint(const double*);

    void copyMeshFrom(const Mesh &);

    void makeIcosahedron();
    void makeSphere(int);
    void makeBiConcave(int);

    void makeOblate(int, double);
    void makeProlate(int, double);

    void buildEdgeList();
    void checkOrientation();

    void calcDoubleLayerLinearTerm(double*);
    void calcDoubleLayerJump();

    // These functions are generic, so include them here
    void tensionForce(const DoubArray1d &, DoubArray2d &);
    void velDiv(const DoubArray2d &, DoubArray1d &);

    void readHDF5(const char *);
    void readDATsingle(const char *);
    void writeHDF5(const char *);
    void writeTecplot(const char*);
};

void Mesh_vertOneRingNbrs(int nvert, int nface, const int (*f2v)[3],
	vector<vector<int> > &nbrs);

void Mesh_edgeList(int nvert, int nface, const int (*f2v)[3], 
	int &nedge, int (*e2v)[2], int (*e2f)[2]);

void Mesh_faceNbrList(int nvert, int nface, const int (*f2v)[3], 
	int nedge, const int (*e2v)[2], const int (*e2f)[2], 
	int (*f2e)[3], int (*f2f)[3]);

void Mesh_reorderF2v(int nvert, int nface, int (*f2v)[3]);

void Mesh_edgeSwap(Mesh&);

void Mesh_slice(const Mesh&, double zmin, double zmax, int nz, MArray<double,1> &vol_slice);

#endif
