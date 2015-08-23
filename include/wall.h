#ifndef WALL_H_
#define WALL_H_

#include "classdec.h"
#include "mesh.h"

class Wall : public Mesh
{
public:
    Wall();
    ~Wall();

    double velTar[3];		// target velocity

    MArray<double,2> f;		// force density
};


#endif 
