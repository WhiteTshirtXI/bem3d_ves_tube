// Mesh edges
#ifndef EDGE_H
#define EDGE_H

#include "classdec.h"

// Note:
//   As one moves from vert[0] to vert[1], face[0] is on the left hand side
struct Edge {
    // ivert -- index to the two end points
    // iface -- index to the two neighboring faces
    int ivert[2], iface[2];
    Point *vert[2];
    Tri *face[2];

    Edge();
    Edge(const Edge &);
    ~Edge();

    double length();
    double bendAngle();
    void gradBendAngle(double, double (*)[3], double (*)[3]);
};

#endif
