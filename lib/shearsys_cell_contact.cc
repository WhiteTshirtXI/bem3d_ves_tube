#include "cxxheaders.h"
#include "shearsys.h"
#include "collision.h"
#include "ewald.h"
#include "mblas.h"
#include "param.h"
#include "geom_oper.h"
#include "debugfunc.h"

/* ShearSys::cellNoContact */


/* Prevent cells from contacting each other
 * Arguments:
 *   DIST_EPS -- distance thresh hold */
void ShearSys::cellNoContact(double DIST_EPS)
{
    vector<Point*> &vlist = vertList[CELL];
    vector<NbrList*> nlists;
    nlists.push_back(&nlist_phys[CELL][CELL]);
    nlists.push_back(&nlist_phys[CELL][RIGID]);

    collision::forceSeparation(vlist, nlists, DIST_EPS);
}


/* Check whether a point lies within a cell 
 * Argument:
 *   x0 -- the coordinate of the point
 *   whichCell -- index of the cell that the point penetrates 
 * Algorithm:
 *   -- If a point is outside a cell, the total spherical angle is 0 */
bool ShearSys::pointInsideSomeCell(const double *x0, int *whichCell)
{
    bool is_interior = false;
    if (whichCell) *whichCell = -1;

    for (int icell = 0; icell < numCells(); icell++) {
	Cell &cell = cells[icell];

	double xtmp[3];
	FOR_I3 xtmp[i] = x0[i];
	ewald::to_CloseBy(cell.center, xtmp);

	double xmin[3], xmax[3];
	cell.getCoordRange(xmin, xmax);

	if (   xtmp[0] < xmin[0] || xtmp[0] > xmax[0] 
	    || xtmp[1] < xmin[1] || xtmp[1] > xmax[1]
	    || xtmp[2] < xmin[2] || xtmp[2] > xmax[2] ) continue;

	double sangle = 0.0;
	for (int iface = 0; iface < cell.numFaces(); iface++) {
	    Tri &face = cell.faces[iface];
	    double xtri[3][3];
	    FOR_I3 FOR_J3 xtri[i][j] = face.vert[i]->x[j];
	    sangle += tri_solidAngle(xtmp, xtri);
	}
	if (fabs(sangle) > 1.E-2) {
	    is_interior = true;
	    if (whichCell) *whichCell = icell;
	    break;
	}
    } // icell

    return is_interior;
}
