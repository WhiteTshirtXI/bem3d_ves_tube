#include "cxxheaders.h"
#include "mblas.h"
#include "point.h"
#include "tri.h"
#include "edge.h"
#include "geom_oper.h"

/* Edge::Edge 
 * Edge::~Edge 
 * Edge::length
 * Edge::bendAngle 
 * Edge::gradBendAngle */

Edge::Edge()
{
    ivert[0] = ivert[1] = -1;
    iface[0] = iface[1] = -1;

    vert[0] = vert[1] = NULL;
    face[0] = face[1] = NULL;
}


Edge::Edge(const Edge &ed)
{
    ivert[0] = ed.ivert[0];
    ivert[1] = ed.ivert[1];

    iface[0] = ed.iface[0];
    iface[1] = ed.iface[1];

    vert[0] = vert[1] = NULL;
    face[0] = face[1] = NULL;
}


Edge::~Edge()
{}


double Edge::length()
{
    double xx[3];
    FOR_I3 xx[i] = vert[1]->x[i] - vert[0]->x[i];
    return sqrt(xx[0]*xx[0] + xx[1]*xx[1] + xx[2]*xx[2]);
}


// See the comments of bendAngle in geom_oper.cc
double Edge::bendAngle()
{
    Point *p0, *p1, *p2, *p3;

    p0 = vert[0];
    p1 = vert[1];

    for (int l = 0; l < 3; ++l) {
        p2 = face[0]->vert[l];
	if (p2 != p0 && p2 != p1) break;
    }

    for (int l = 0; l < 3; ++l) {
        p3 = face[1]->vert[l];
	if (p3 != p0 && p3 != p1) break;
    }

    return ::bendAngle(p0->x, p1->x, p2->x, p3->x);
}


/* The gradient of bending angle w.r.t. vertex coordinates
 * Arguments:
 *   th -- the bending angle (input)
 *   D0,D1 -- the gradient tensor
 * Note:
 *  -- D[i][j] = dth/dx[i][j] for vertices in face[0] and face[1]
 *  -- For the two vertices shared by the two neighboring faces, the gradients 
 *     are the sum of the corresponding components in D0 and D1 */
void Edge::gradBendAngle(double th, double (*D0)[3], double (*D1)[3])
{
    double *n0, *n1;
    double dn0[3][3][3], dn1[3][3][3];
    
    // Get the normal of the two faces and their gradients wrt vertex coordinates
    n0 = face[0]->normal;
    n1 = face[1]->normal;

    face[0]->gradNormal(dn0);
    face[1]->gradNormal(dn1);

    if (fabs(th) > 0.25*M_PI) {
        // non-coplanar case
	// dcosth = -sinth * dth
	
	// costh = n0 * n1
	// dcosth = dn0 * n1 + n0 * dn1
	for (int i = 0; i < 3; ++i)
	for (int ii = 0; ii < 3; ++ii) {
	    D0[i][ii] = 0.0;
	    D1[i][ii] = 0.0;

	    for (int jj = 0; jj < 3; ++jj) {
		D0[i][ii] += dn0[jj][i][ii] * n1[jj];
		D1[i][ii] += dn1[jj][i][ii] * n0[jj];
	    }
	}

	// converg d(costh)/dx to d(th)/dx
	double fac = -1.0/sin(th);
	for (int i = 0; i < 3; ++i)
	for (int ii = 0; ii < 3; ++ii) {
	    D0[i][ii] *= fac;
	    D1[i][ii] *= fac;
	}
    }
    else {
        // when the two faces are almost co-planar
	// dsinth = costh * dth
	
	// sinth = (n0 x n1) * t
	// dsinth = (dn0 x n1)*t + (n0 x dn1)*t 
	//        = dn0 * (n1 x t) + dn1 * (t x n0)
        double t[3], tt;
	double n1xt[3], txn0[3];

	for (int ii = 0; ii < 3; ++ii) 
	    t[ii] = vert[1]->x[ii] - vert[0]->x[ii];
        tt = m_dnrm2(3, t);
	m_dscal(3, 1.0/tt, t);

	cross_product(n1, t, n1xt);
	cross_product(t, n0, txn0);

	for (int i = 0; i < 3; ++i)
	for (int ii = 0; ii < 3; ++ii) {
	    D0[i][ii] = 0.0;
	    D1[i][ii] = 0.0;

	    for (int kk = 0; kk < 3; ++kk) {
		D0[i][ii] += dn0[kk][i][ii] * n1xt[kk];
		D1[i][ii] += dn1[kk][i][ii] * txn0[kk];
	    }
	}

	// change d(sinth)/dx to d(th)/dx
	double fac = 1.0/cos(th);
	for (int i = 0; i < 3; ++i)
	for (int ii = 0; ii < 3; ++ii) {
	    D0[i][ii] *= fac;
	    D1[i][ii] *= fac;
	}
    }
}
