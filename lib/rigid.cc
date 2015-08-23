#include "cxxheaders.h"
#include "rigid.h"
#include "mathfunc.h"
#include "geom_oper.h"
#include "quadrature.h"
#include "ewald.h"

/* Rigid::Rigid 
 * Rigid::~Rigid 
 * void Rigid::projectRigidVel
 * void Rigid::move */

Rigid::Rigid()
{ }


Rigid::~Rigid()
{ }


/* Project a velocity vector to rigid body motion
 * Argument:
    u -- the velocity to be projected */
void Rigid::projectRigidVel(MArray<double,2> &u)
{
    // calculate the translating and rotating velocity by projection
    double UT[3], OMG[3];
    calcTransRotatVel(u, UT, OMG);

    // Do the projection
    for (int ivert = 0; ivert < numVerts(); ivert++) {
	double xx[3], omgxx[3];

	FOR_I3 xx[i] = verts[ivert].x[i] - center[i];
	cross_product(OMG, xx, omgxx);

	// Add translation
	FOR_I3 u(ivert,i) = UT[i] + omgxx[i];
    } // ivert
}


/* Move the rigid body, assuming the velocity is projected from the _g array
 * Argument:
 *   u -- the velocity
 *   Ts -- the time step */
void Rigid::move(const MArray<double,2> &u, double dt)
{
    double UT[3], dx[3];
    double omg[3], rot_axis[3], rot_angle;

    // Calculate translating and rotation velocity
    calcTransRotatVel(u, UT, omg);

    // Translating distance
    FOR_I3 dx[i] = dt*UT[i];

    // Rotating axis and rotating angle
    double omg2 = m_dnrm2(3, omg);

    if (omg2 > 1.E-10) {
        FOR_I3 rot_axis[i] = omg[i]/omg2;
	rot_angle = dt*omg2;
    } else {
        rot_axis[0] = 1.0;
        rot_axis[1] = 0.0;
        rot_axis[2] = 0.0;
	rot_angle = 0.0;
    }

    for (int ivert = 0; ivert < numVerts(); ivert++) {
        Point &vert = verts[ivert];

	double xx[3];
	FOR_I3 xx[i] = vert.x[i] - center[i];
	rotate_vector(rot_axis, rot_angle, xx);

	FOR_I3 vert.x[i] = center[i] + xx[i] + dx[i];
    } // ivert
}
