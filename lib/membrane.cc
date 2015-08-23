#include "cxxheaders.h"
#include "membrane.h"
#include "mesh.h"
#include "mathfunc.h"

/* void tri_strain_energy 
 * void edge_bend_energy  */


/* Compute the strain energy and its derivatives on a triangle 
 * Arguments: 
 *  T -- current triangle 
 *  Tref -- reference triangle 
 *  ES -- shear modulus 
 *  ED -- dilatation modulus 
 *  W -- the strain energy 
 *  DW[i][j] = DW/Dx[i][j] 
 * Note: 
 *  -- W = Es/2*(I1 - I2 + I1^2/4) + ED/8*I2*I2 
 *  -- shear modulus = ES
 *  -- dilatational modulus = ED */
void tri_strain_energy(Tri &T, Tri &Tref, double ES, double ED,
		double *W, double (*DW)[3])
{
    // tangents
    double a0[3], a1[3];
    for (int ii = 0; ii < 3; ++ii) {
	a0[ii] = T.vert[0]->x[ii] - T.vert[2]->x[ii];
	a1[ii] = T.vert[1]->x[ii] - T.vert[2]->x[ii];
    }

    double I1 = m_ddot(4, &T.a[0][0], &Tref.arcp[0][0]) - 2; 
    double I2 = T.detA/Tref.detA - 1;

    if (W) {
	*W = 0.5*ES*(I1 - I2 + 0.25*I1*I1) + 0.125*ED*I2*I2;
	*W *= Tref.area;
    }

    if (DW) {
	double W1, W2;
	double c;

	W1 = 0.5*ES*(1 + 0.5*I1);
	W2 = -0.5*ES + 0.25*ED*I2;

	for (int i = 0; i < 2; ++i) {
	    fill_n(DW[i], 3, 0.0);

	    c = W1*2*Tref.arcp[i][0] + W2*2*(T.detA/Tref.detA)*T.arcp[i][0];
	    m_daxpy(3, c, a0, DW[i]);

	    c = W1*2*Tref.arcp[i][1] + W2*2*(T.detA/Tref.detA)*T.arcp[i][1];
	    m_daxpy(3, c, a1, DW[i]);

	    // Scale by triangle area
	    m_dscal(3, Tref.area, DW[i]);
	} // i

	for (int ii = 0; ii < 3; ++ii) {
	    DW[2][ii] = -DW[0][ii] - DW[1][ii];
	}
    }
}


/* Calculate the bending energy at an edge 
 * Arguments: 
 *  ed -- edge
 *  sa -- spontaneous angle
 *  c -- edge bending coefficient
 *  W -- bending energy  
 *  DW0, DW1 -- derivatives of bending energy w.r.t the 
 *  		1st and 2nd triangle's coordinates */
void edge_bend_energy(Edge &ed, double sa, double c, 
		double *W, double (*DW0)[3], double (*DW1)[3])
{
    double th = ed.bendAngle();;

    if (W) {
        *W = c*(1 - cos(th - sa));
    }

    if (DW0 && DW1) {
        ed.gradBendAngle(th, DW0, DW1);

	double fac = c*sin(th - sa);
	m_dscal(9, fac, *DW0);
	m_dscal(9, fac, *DW1);
    }
}
