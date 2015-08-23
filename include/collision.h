#ifndef COLISION_H
#define COLISION_H

#include "ewald.h"

namespace collision
{
    // Collision information
    struct Cld {
        int ivert;
        double spr, dx[3];

        Cld() : ivert(-1), spr(1.E10)
        {
            dx[0] = dx[1] = dx[2] = 0.0;
        }
    };

    // Compare two collision
    inline bool CldCmp(const Cld &c1, const Cld &c2) {
        return ( (c1.ivert < c2.ivert) || 
                 (c1.ivert == c2.ivert && c1.spr < c2.spr) );
    }

    void updateCldInfo(NbrList &nlist, double DIST_EPS, 
                MArray<int,1> &pcld, vector<Cld> &clds);

    void forceSeparation(vector<Point*> &vlist, vector<NbrList*> &nlists, double DIST_EPS);
}

#endif
