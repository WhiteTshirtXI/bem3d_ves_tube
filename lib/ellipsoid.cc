#include "cxxheaders.h"
#include "ellipsoid.h"
#include "mathfunc.h"

/* void ellipsoid::calcLatLonAngles
 * void ellipsoid::calcLocalGeom */


void ellipsoid::calcLatLonAngles(double A, double B, double C, const double *x,
		double &th, double &phi)
{
    assert(fabs(x[2]) <= C);

    th = acos(x[2]/C);
    phi = atan2(x[1]/B, x[0]/A);
}


void ellipsoid::calcLocalGeom(double A, double B, double C, double th, double phi,
    		double *x, double *nrml, double &H, double &K)
{
    x[0] = A*sin(th)*cos(phi);
    x[1] = B*sin(th)*sin(phi);
    x[2] = C*cos(th);

    if (fabs(sin(th)) < 1.E-5) {
        FOR_I3 nrml[i] = x[i];
        normalizeVec3D(nrml);

        H = -0.5*(C/(A*A) + C/(B*B));
        K = (C*C)/(A*A*B*B);
    } else {
	double a1[3] = { A*cos(th)*cos(phi), B*cos(th)*sin(phi), -C*sin(th) };
	double a2[3] = { -A*sin(th)*sin(phi), B*sin(th)*cos(phi), 0.0 };

	double a11[3] = { -A*sin(th)*cos(phi), -B*sin(th)*sin(phi), -C*cos(th) };
	double a12[3] = { -A*cos(th)*sin(phi), B*cos(th)*cos(phi), 0.0 };
	double a22[3] = { -A*sin(th)*cos(phi), -B*sin(th)*sin(phi), 0.0 };

        cross_product(a1, a2, nrml);
        normalizeVec3D(nrml);

        double a[2][2], arcp[2][2], b[2][2];
        a[0][0] = m_ddot(3, a1, a1);
        a[0][1] = m_ddot(3, a1, a2);
        a[1][0] = a[0][1];
        a[1][1] = m_ddot(3, a2, a2);

        invMat2(a, arcp);
        invMat2(a, arcp);

        b[0][0] = m_ddot(3, a11, nrml);
        b[0][1] = m_ddot(3, a12, nrml);
        b[1][0] = b[0][1];
        b[1][1] = m_ddot(3, a22, nrml);

        H = 0.5*(arcp[0][0]*b[0][0] + arcp[0][1]*b[0][1]
               + arcp[1][0]*b[1][0] + arcp[1][1]*b[1][1]);

        K = b[0][0]*b[1][1] - b[0][1]*b[1][0];
        K *= (arcp[0][0]*arcp[1][1] - arcp[0][1]*arcp[1][0]);
    }
}
