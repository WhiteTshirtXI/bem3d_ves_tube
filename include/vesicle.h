#ifndef VESICLE_H
#define VESICLE_H

#include "classdec.h"
#include "mesh.h"
#include "subdiv.h"

class Vesicle : public Mesh
{
public:
    double areaTar, volTar;	// target surface area and volume

    MArray<double,1> sigma;	// surface tension
    MArray<double,2> v;		// velocity

    MArray<double,2> fbend;	// bending force
    MArray<double,2> ftens;	// tension force
    MArray<double,2> f;	// fbend+ftens, added by spann for compatibility with blood code
    

    double EB;			// bending modulus
    double H0;			// spontaneous curvature

    Vesicle();
    ~Vesicle();

    vector<int> vertValence;		// vertex valences
    vector<vector<int> > ctrlPts;	// control points

    void buildControlPoints();

	void vesaddtension(); //Added by spann, sets f=fbend+ftens
    double bendEnergy();
    void bendForce(MArray<double,2> &f);
    void diffBendForce(const MArray<double,2> &dx, MArray<double,2> &df);

    void calcMeshRelaxVel(MArray<double,2> &v);
    int relaxMesh(int niter = 1, double strength = 1.E-2);

    void gradVolArea(MArray<double,1> &gradV, MArray<double,1> &gradA);
    void adjustVolArea(double Vtar, double Atar);

    //**********************************************************************
    // Quadrature points
    struct QuadPoint {
        double s, t, wght;

        double x[3], dx[2][3], ddx[2][2][3];
        double a[2][2], b[2][2];
        double dx_rcp[2][3], a_rcp[2][2], b_rcp[2][2];
	double nrml[3], area, H, K;
    };

    void calcQuadPoint(
    	const DoubArray1d &f, const DoubArray2d &df, const DoubArray3d &ddf, 
	const DoubArray2d &xctrl, QuadPoint &q);

    //**********************************************************************
};

#endif
