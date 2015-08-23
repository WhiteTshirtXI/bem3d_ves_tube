#ifndef CELL_H_
#define CELL_H_

#include "classdec.h"
#include "mesh.h"

class Cell : public Mesh
{
public:
    Cell *cellRef;		// Reference cell

    double volTar;		// Target volume

    // ES -- shear modulus
    // ED -- dilatational modulus
    // EB -- bending modulus
    double ES, ED, EB;

    MArray<double,2> felas;	// elastic force
    MArray<double,2> fbend;	// surface force
    MArray<double,2> f;		// total surface force
    MArray<double,2> v;		// velcoity

    Cell();
    ~Cell();

    // Spontaneous bending angle
    double spontaneousAngle() { 
        return sqrt( 16*M_PI/(3*sqrt(3.0)*numFaces()) );
    }

    double elasticEnergy();
    void elasticForce(MArray<double,2> &);

    double bendingEnergy();
    void bendingForce(MArray<double,2> &);

    void calcSurfForce(MArray<double,2> &);
    void diffForce(const MArray<double,2> &, MArray<double,2> &);
};

#endif
