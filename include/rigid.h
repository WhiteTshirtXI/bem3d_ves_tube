#ifndef RIGID_H
#define RIGID_H

#include "classdec.h"
#include "mesh.h"

class Rigid : public Mesh
{
public:
    Rigid();
    ~Rigid();

    MArray<double,2> f;		// surface force
    MArray<double,2> v;		// surface velocity

    void projectRigidVel(MArray<double,2> &);
    void move(const MArray<double,2> &, double);
};

#endif
