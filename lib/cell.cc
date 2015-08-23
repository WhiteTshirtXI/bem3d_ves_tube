#include "cxxheaders.h"
#include "cell.h"
#include "edge.h"
#include "membrane.h"
#include "mathfunc.h"

/* Cell::Cell 
 * Cell::~Cell 
 * Cell::calcSurfForce 
 * Cell::elasticEnergy 
 * Cell::elasticForce 
 * Cell::bendingEnergy 
 * Cell::bendingForce 
 * Cell::diffForce */

Cell::Cell()
{ 
    cellRef = NULL;
}


Cell::~Cell()
{ 
}


void Cell::calcSurfForce(MArray<double,2> &f)
{
    int nvert = numVerts();

    felas.resize(nvert,3);
    fbend.resize(nvert,3);
    f.resize(nvert,3);

    elasticForce(felas);
    bendingForce(fbend);

    f = 0.0;
    f += felas;
    f += fbend;
}


double Cell::elasticEnergy()
{
    double E = 0.0;

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &T = faces[iface];
        Tri &T_ref = cellRef->faces[iface];
        double dE;

        tri_strain_energy(T, T_ref, ES, ED, &dE, NULL);
	E += dE;
    }

    return E;
}


void Cell::elasticForce(MArray<double,2> &f)
{
    int nvert = numVerts();

    // Init
    f.resize(nvert,3);
    f = 0.0;

    for (int iface = 0; iface < numFaces(); ++iface) {
	Tri &T = faces[iface];
	Tri &T_ref = cellRef->faces[iface];

	double DW[3][3];
	tri_strain_energy(T, T_ref, ES, ED, NULL, DW);

	for (int l = 0; l < 3; ++l) {
	    int ivert = T.ivert[l];
	    m_daxpy(3, 1.0/vertArea(ivert), DW[l], &f(ivert,0));
	}
    } // iface
}


double Cell::bendingEnergy()
{
    // Init
    double E = 0.0;

    const double thRef = spontaneousAngle();
    const double kappa = 2.0*sqrt(3.0)*EB;

    for (int ied = 0; ied < numEdges(); ++ied) {
        Edge &ed = edges[ied];

	double dE;
	edge_bend_energy(ed, thRef, kappa, &dE, NULL, NULL);

	E += dE;
    } // ied

    return E;
}


void Cell::bendingForce(MArray<double,2> &f)
{
    int nvert = numVerts();
    int nface = numFaces();

    // Init
    f.resize(nvert,3);
    f = 0.0;

    const double thRef = spontaneousAngle();
    const double kappa = 2*sqrt(3.0)*EB;

    for (int ied = 0; ied < numEdges(); ied++) {
	Edge &ed = edges[ied];

	double DW0[3][3], DW1[3][3];

	edge_bend_energy(ed, thRef, kappa, NULL, DW0, DW1);

	for (int l = 0; l < 3; l++) {
	    int ivert = ed.face[0]->ivert[l];
	    m_daxpy(3, 1.0/vertArea(ivert), DW0[l], &f(ivert,0));

	    ivert = ed.face[1]->ivert[l];
	    m_daxpy(3, 1.0/vertArea(ivert), DW1[l], &f(ivert,0));
	}
    } // ied
}


/* The change to surface force due to surface perturbation
 * Arguments:
 *  dx -- surface perturbation
 *  df -- change to bending force 
 * Note:
 *   -- Assume the original f(:,:) is already calcualted */
void Cell::diffForce(const MArray<double,2> &dx, MArray<double,2> &df)
{
    const int nvert = numVerts();
    const int nface = numFaces();

    // Save old coordinates
    MArray<double,2> xOld(nvert,3);
    getCoords(xOld);

    // Determine the scaling 
    const double meshSize = sqrt( 4*area/(nface*sqrt(3.0)) );

    double max_dx = 0.0;
    for (int ivert = 0; ivert < nvert; ivert++) {
        FOR_J3 max_dx = max(max_dx, fabs(dx(ivert,j)) );
    }

    // Trivial case
    if (max_dx < 1.E-10*meshSize) {
        df.resize(nvert,3);
	df = 0.0;
	return;
    }

    // Perturb the surface and calculate the new force
    const double scal = 1.E-3*meshSize/max_dx;
    for (int ivert = 0; ivert < nvert; ivert++) {
        Point &vert = verts[ivert];
	m_daxpy(3, scal, &dx(ivert,0), vert.x);
    }
    updateGeometry();

    MArray<double,2> felas_new(nvert,3);
    MArray<double,2> fbend_new(nvert,3);

    elasticForce(felas_new);
    bendingForce(fbend_new);

    // Restore the original state
    setCoords(xOld);
    updateGeometry();

    // Use finite difference
    df.resize(nvert,3);

    df = 0.0;
    df += felas_new;
    df += fbend_new;
    df -= felas;
    df -= fbend;
    df *= 1./scal;
}
