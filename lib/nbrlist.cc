#include "cxxheaders.h"
#include "point.h"
#include "tri.h"
#include "ewald.h"
#include "mathfunc.h"
#include <set>
using std::set;

/* NbrList::build
 * NbrList::build_rect
 * NbrList::build_general
 * NbrList::findNumNbrPoints */

using ewald::phys::rc;

/* Build neighbor list
 * Argument:
 *   points -- the target point list
 *   slist -- the source list 
 *   NO_SELF -- if true, do not include self-interactions (i.e. between points 
 *   		and faces on the same mesh)
 *   ONLY_SELF -- if true, only include self-interactions */
void NbrList::build(const vector<Point*> &points, SrcList &slist,
                bool NO_SELF, bool ONLY_SELF)
{
    if (ewald::bctype == ewald::RECT)
        build_rect(points, slist, NO_SELF, ONLY_SELF);
    else
        build_general(points, slist, NO_SELF, ONLY_SELF);
}



void NbrList::build_rect(const vector<Point*> &points, SrcList &slist,
                bool NO_SELF, bool ONLY_SELF)
{
    // Init
    verts.clear();
    faces.clear();
    dists.clear();
    firstNbr.clear();

    if (points.size() == 0 || slist.numFaces() == 0) return;

    for (int i = 0; i < points.size(); i++) {
        Point *vert = points[i];
        int i0, i1, i2;
        slist.getBlkNum(vert->x, i0, i1, i2);

        if (i0 < -1 || i0 > slist.nblk[0] ||
            i1 < -1 || i1 > slist.nblk[1] ||
            i2 < -1 || i2 > slist.nblk[2] ) continue;

        verts.push_back(vert);
        firstNbr.push_back(faces.size());

        for (int j0 = max(i0-1,-1); j0 <= min(i0+1,slist.nblk[0]); j0++)
        for (int j1 = max(i1-1,-1); j1 <= min(i1+1,slist.nblk[1]); j1++)
        for (int j2 = max(i2-1,-1); j2 <= min(i2+1,slist.nblk[2]); j2++) {
            for (int j = slist.hoc(j0,j1,j2); j >= 0; j = slist.next(j)) {
                Tri *tri = slist.faces[j];

                // Exclude points on the same mesh
                if (NO_SELF) {
                    if (vert->mesh != NULL && vert->mesh == tri->mesh) continue;
                }

                if (ONLY_SELF) {
                    if (vert->mesh != NULL && vert->mesh != tri->mesh) continue;
                }

                double xtar[3];
                FOR_K3 xtar[k] = vert->x[k];
                ewald::to_CloseBy(tri->xc, xtar);
                double rr = tri->minDistToPoint_Approx(xtar);

                if (rr < rc) {
                    faces.push_back(tri);
                    dists.push_back(rr);
                }
            } // j
        } // j0,j1,j2
    } // i

    firstNbr.push_back(faces.size());
}


void NbrList::build_general(const vector<Point*> &points, SrcList &slist,
                bool NO_SELF, bool ONLY_SELF)
{
    using ewald::L;
    using ewald::iL;

    // Init
    verts.clear();
    faces.clear();
    dists.clear();
    firstNbr.clear();

    if (points.size() == 0 || slist.numFaces() == 0) return;

    double xminBuf[3], xmaxBuf[3];
    FOR_D3 {
        xminBuf[d] = slist.lb[d] - rc;
	xmaxBuf[d] = slist.ub[d] + rc;
    }

    // Find the bounding box in the lattice coordinates
    double sminBuf[3], smaxBuf[3];
    FOR_D3 sminBuf[d] = FLT_MAX;
    FOR_D3 smaxBuf[d] = -FLT_MAX;

    for (int ix = 0; ix <= 1; ix++)
    for (int iy = 0; iy <= 1; iy++)
    for (int iz = 0; iz <= 1; iz++) {
        double x[3];
	x[0] = (ix == 0) ? xminBuf[0] : xmaxBuf[0];
	x[1] = (iy == 0) ? xminBuf[1] : xmaxBuf[1];
	x[2] = (iz == 0) ? xminBuf[2] : xmaxBuf[2];

	double r[3];
	ewald::phys_to_lattice_coord(x, r);

	FOR_D3 {
	    sminBuf[d] = min( sminBuf[d], r[d] );
	    smaxBuf[d] = max( smaxBuf[d], r[d] );
        }
    }

    for (int i = 0; i < points.size(); i++) {
        Point *vert = points[i];

	// Find all the image points that lies in the extended
	// bounding box
	vector<double> ximgs;

	double sbase[3];	// base-line lattice coordinate
	ewald::phys_to_lattice_coord( vert->x, sbase );
	FOR_D3 sbase[d] += ceil(sminBuf[d] - sbase[d]);

	double s[3];
	for ( s[0] = sbase[0]; s[0] < smaxBuf[0]; s[0]++ )
	for ( s[1] = sbase[1]; s[1] < smaxBuf[1]; s[1]++ )
	for ( s[2] = sbase[2]; s[2] < smaxBuf[2]; s[2]++ ) {
	    double x[3];
	    ewald::lattice_to_phys_coord(s, x);

	    if ( x[0] > xminBuf[0] && x[0] < xmaxBuf[0] &&
	         x[1] > xminBuf[1] && x[1] < xmaxBuf[1] &&
		 x[2] > xminBuf[2] && x[2] < xmaxBuf[2] ) {
		ximgs.push_back(x[0]);
		ximgs.push_back(x[1]);
		ximgs.push_back(x[2]);
	     }
	}

	if (ximgs.size() == 0) continue;

	verts.push_back(vert);
	firstNbr.push_back(faces.size());

	for (int iimg = 0; iimg < ximgs.size()/3; iimg++) {
	    double xtar[3];
	    xtar[0] = ximgs[iimg*3];
	    xtar[1] = ximgs[iimg*3+1];
	    xtar[2] = ximgs[iimg*3+2];

	    int i0, i1, i2;
	    slist.getBlkNum(xtar, i0, i1, i2);

	    if (i0 < -1 || i0 > slist.nblk[0] ||
		i1 < -1 || i1 > slist.nblk[1] ||
		i2 < -1 || i2 > slist.nblk[2] ) continue;

	    for (int j0 = max(i0-1,0); j0 <= min(i0+1,slist.nblk[0]-1); j0++)
	    for (int j1 = max(i1-1,0); j1 <= min(i1+1,slist.nblk[1]-1); j1++)
	    for (int j2 = max(i2-1,0); j2 <= min(i2+1,slist.nblk[2]-1); j2++) {
		for (int j = slist.hoc(j0,j1,j2); j >= 0; j = slist.next(j)) {
		    Tri *tri = slist.faces[j];

		    // Exclude points on the same mesh
		    if (NO_SELF) {
			if (vert->mesh != NULL && vert->mesh == tri->mesh) continue;
		    }

		    if (ONLY_SELF) {
			if (vert->mesh != NULL && vert->mesh != tri->mesh) continue;
		    }

		    double rr = tri->minDistToPoint_Approx(xtar);

		    if (rr < rc) {
			faces.push_back(tri);
			dists.push_back(rr);
		    }
		} // j
	    } // j0,j1,j2
	}
    } // i

    firstNbr.push_back(faces.size());
}


/* Build neighbor list without hashing */
void NbrList::build(const vector<Point*> &tgts, const vector<Tri*> &srcs)
{
    // Init
    verts = tgts;
    faces.clear();
    dists.clear();
    firstNbr.clear();

    for (int i = 0; i < verts.size(); i++) {
        firstNbr.push_back(faces.size());

        for (int j = 0; j < srcs.size(); j++) {
            Tri *tri = srcs[j];

            double xtar[3];
            FOR_K3 xtar[k] = verts[i]->x[k];
            ewald::to_CloseBy(tri->xc, xtar);
            double rr = tri->minDistToPoint_Approx(xtar);

            if (rr < rc) {
                faces.push_back(tri);
                dists.push_back(rr);
            }
        } // j
    } // i

    firstNbr.push_back(faces.size());
}


/* Find the interacting vertices
 * Argument:
 *  nnbr -- the number of interacting source points of each target point */
void NbrList::findNumNbrPoints(int *nnbr)
{
    set<int> nbr_index;

    for (int i = 0; i < verts.size(); i++) {
        nbr_index.clear();
        for (int p = firstNbr[i]; p < firstNbr[i+1]; p++) {
            Tri *tri = faces[p];
            for (int l = 0; l < 3; l++)
                nbr_index.insert(tri->vert[l]->Gindx);
        }
        nnbr[i] = nbr_index.size();
    } // i
}
