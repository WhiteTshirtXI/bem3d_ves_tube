#ifndef REDCELL_H
#define REDCELL_H

#include "classdec.h"
#include "mesh.h"

class RedCell : public Mesh
{
public:
    RedCell *cellRef;		// reference cell
    double volTar;		// target volume
    double areaTar;		// target area

    double ES;		// shear modulus
    double ED;		// dilatational modulus
    double EB;		// bending modulus
    double H0;		// local spontaneous curvature
    double EB_GLB;	// global bending modulus
    double H0_GLB;	// global spontaneous curvature
    double muM;		// membrane viscosity

    double sigma;	// surface tension
    // Error in surface area
    double area_E0, area_E1, area_E2;
    // PID controller parameters for surface tension
    static const double KP = 1.0;
    static const double KI = 5.0;
    static const double KD = 1.0;

    int elas_model; 
    enum { SKALAK_MODEL, PNAS2002_MODEL };

    // Force and velcoity
    MArray<double,2> f;	
    MArray<double,2> v;		

    RedCell();
    ~RedCell();

    virtual void updateGeometry();

    double tri_elasticity(Tri &T, Tri &Tref, double (*grad)[3]=NULL);
    double tri_elasticity_Skalak(Tri &T, Tri &Tref, double (*grad)[3]=NULL);
    double tri_elasticity_PNAS2002(Tri &T, Tri &Tref, double (*grad)[3]=NULL);

    double elasticEnergy();
    void elasticForce(MArray<double,2> &);

    double bendEnergy();
    void bendForce(MArray<double,2> &);

    void diffElasBendForce(MArray<double,2> &, MArray<double,2> &);

    void tensionForce(double, MArray<double,2> &);

    double viscousForce(MArray<double,2> &, MArray<double,2> &);
    void strainRate(MArray<double,2> &, MArray<double,1> &);

    void tweezerForce(double, MArray<double,2> &);

    //**********************************************************************
    // Use Loop subdivision for bending-force calculation
    vector<int> vertValence;		// vertex valences
    vector<vector<int> > ctrlPts;	// control points
    void buildControlPoints();

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
