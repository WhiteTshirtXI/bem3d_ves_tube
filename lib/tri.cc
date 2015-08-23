#include "cxxheaders.h"
#include "point.h"
#include "tri.h"
#include "mathfunc.h"

/* Tri::Tri 
 * ~Tri::Tri 
 * Tri::updateGeometry 
 * Tri::baryCoord
 * Tri::minDistToPoint_Approx 
 * Tri::gradNormal 
 * Tri::gradArea 
 * Tri::secondDerivArea 
 * Tri::minAngle */


// default constructor
Tri::Tri() 
{
    ivert[0] = ivert[1] = ivert[2] = -1;
    vert[0] = vert[1] = vert[2] = NULL;
    mesh = NULL;
}


Tri::Tri(int i0, int i1, int i2) 
{
    ivert[0] = i0;
    ivert[1] = i1;
    ivert[2] = i2;

    vert[0] = vert[1] = vert[2] = NULL;
    mesh = NULL;
}


// copy constructor
Tri::Tri(const Tri& t) 
{
    ivert[0] = t.ivert[0];
    ivert[1] = t.ivert[1];
    ivert[2] = t.ivert[2];

    vert[0] = vert[1] = vert[2] = NULL;
    mesh = NULL;
}


// destructor
Tri::~Tri() 
{ }


void Tri::updateGeometry() 
{
    double a0[3], a1[3], idetA;

    for (int ii = 0; ii < 3; ++ii) {
	// centroid
        xc[ii] = THRD*(vert[0]->x[ii] + vert[1]->x[ii] + vert[2]->x[ii]);

	// tangent
        a0[ii] = vert[0]->x[ii] - vert[2]->x[ii];
        a1[ii] = vert[1]->x[ii] - vert[2]->x[ii];
    }

    // Metric tensor and its inverse
    a[0][0] = m_ddot(3, a0, a0);
    a[0][1] = a[1][0] = m_ddot(3, a0, a1);
    a[1][1] = m_ddot(3, a1, a1);

    detA = a[0][0]*a[1][1] - a[0][1]*a[1][0];
    idetA = 1.0/detA;

    arcp[0][0] = idetA*a[1][1];
    arcp[1][1] = idetA*a[0][0];
    arcp[0][1] = arcp[1][0] = -idetA*a[0][1];

    // Surface normal and area
    normal[0] = a0[1]*a1[2] - a0[2]*a1[1];
    normal[1] = a0[2]*a1[0] - a0[0]*a1[2];
    normal[2] = a0[0]*a1[1] - a0[1]*a1[0];

    detJ = m_dnrm2(3, normal);
    m_dscal(3, 1.0/detJ, normal);

    area = 0.5*detJ;

    // Calculate _normal and _offset arrays
    // These are used for calculate the approximate distance 
    // between the triangle and a point
    m_dcopy(3, normal, _normal[3]);
    _offset[3] = m_ddot(3, _normal[3], xc);

    for (int i = 0; i < 3; i++) {
        int i1 = (i + 1)%3;
        double xx[3];
	FOR_K3 xx[k] = vert[i1]->x[k] - vert[i]->x[k];

	cross_product(xx, normal, _normal[i]);
	normalizeVec3D(_normal[i]);
	_offset[i] = m_ddot(3, _normal[i], vert[i]->x);
    }

    // Calculate the circumcircle's diameter
    double L0 = sqrt(a[0][0]);
    double L1 = sqrt(a[1][1]);
    double L2 = sqrt(fabs(a[0][0] + a[1][1] - 2*a[0][1]));
    diam = max(L0, max(L1, L2));
}


/* Calculate the Barycentric coordinate
 * Arguments:
 *   x -- the coordinate in physical space
 *   s, t -- the coordinate in the reference triangle
 * Note:
 *   1. The reference triangle is (0,0), (1,0), (0,1)
 *   2. Barycentric coordinates are w0=1-s-t, w1=s, w2=t
 */
void Tri::baryCoord(const double *x, double &s, double &t) const
{
    double a1[3], a2[3], dx[3];
    double lhs[2][2], rhs[2];
    double ilhs[2][2];

    FOR_I3 a1[i] = vert[1]->x[i] - vert[0]->x[i];
    FOR_I3 a2[i] = vert[2]->x[i] - vert[0]->x[i];
    FOR_I3 dx[i] = x[i] - vert[0]->x[i];

    lhs[0][0] = m_ddot(3, a1, a1);
    lhs[0][1] = m_ddot(3, a1, a2);
    lhs[1][0] = lhs[0][1];
    lhs[1][1] = m_ddot(3, a2, a2);

    rhs[0] = m_ddot(3, a1, dx);
    rhs[1] = m_ddot(3, a2, dx);

    invMat2(lhs, ilhs);

    s = ilhs[0][0]*rhs[0] + ilhs[0][1]*rhs[1];
    t = ilhs[1][0]*rhs[0] + ilhs[1][1]*rhs[1];
}


// Approximate distance to a point
// Argument:
//  x -- the target point
// Note:
//  This is a conservative estimation, i.e. the approximate distance
//  is smaller than the exact value.
double Tri::minDistToPoint_Approx(const double *x) const
{
    double r0, r1, r2, r3;

    r0 = m_ddot(3, _normal[0], x) - _offset[0];
    r1 = m_ddot(3, _normal[1], x) - _offset[1];
    r2 = m_ddot(3, _normal[2], x) - _offset[2];
    r3 = abs(m_ddot(3, _normal[3], x) - _offset[3]);

    return max(max(r0,r1), max(r2,r3));
}


// Calculate the gradient of the surface normal WRT vertex coordinates
// Arguments:
//   dn[i][j][k] = Dn(i)/Dx(j,k)
void Tri::gradNormal(double (*dn)[3][3]) const
{
    double a0[3], a1[3], a0rcp[3], a1rcp[3];

    for (int ii = 0; ii < 3; ++ii) {
        a0[ii] = vert[0]->x[ii] - vert[2]->x[ii];
        a1[ii] = vert[1]->x[ii] - vert[2]->x[ii];
    }

    for (int ii = 0; ii < 3; ++ii) {
	a0rcp[ii] = arcp[0][0]*a0[ii] + arcp[0][1]*a1[ii];
	a1rcp[ii] = arcp[1][0]*a0[ii] + arcp[1][1]*a1[ii];
    }

    for (int ii = 0; ii < 3; ++ii)
    for (int jj = 0; jj < 3; ++jj) {
        double v0 = -a0rcp[ii]*normal[jj];
	double v1 = -a1rcp[ii]*normal[jj];

        dn[ii][0][jj] = v0;
	dn[ii][1][jj] = v1;
	dn[ii][2][jj] = -(v0 + v1);
    }
}


/* Calculate the gradient of area w.r.t. vertex coordinates 
 * Arguments:
 *   dA[j][k] = dA/dx(j,k) */
void Tri::gradArea(double (*dA)[3]) const
{
    double a0[3], a1[3], a0rcp[3], a1rcp[3];

    FOR_I3 {
        a0[i] = vert[0]->x[i] - vert[2]->x[i];
        a1[i] = vert[1]->x[i] - vert[2]->x[i];
    }

    FOR_I3 {
	a0rcp[i] = arcp[0][0]*a0[i] + arcp[0][1]*a1[i];
	a1rcp[i] = arcp[1][0]*a0[i] + arcp[1][1]*a1[i];
    }

    FOR_I3 {
	dA[0][i] = area*a0rcp[i];
	dA[1][i] = area*a1rcp[i];
	dA[2][i] = -dA[0][i] - dA[1][i];
    }
}


/* Calculate the second derivative of area w.r.t. vertex coordinates
 * Arguments:
 *   D2[p][i][q][j] = d^2 A/dx(p,i) dx(q,j) */
void Tri::secondDerivArea(double (*D2)[3][3][3]) const
{
    // Calculate first derivative 
    double D1[3][3];
    gradArea(D1);

    // Initialize
    m_dclear(81, D2[0][0][0]);

    // Calculate the upper left of the matrix
    for (int l = 0; l < 3; l++)
    for (int m = 0; m < 3; m++) {
        double D2tmp[3][3] = { 0.0 };

        if (l > m) continue;

	if (l == m) {
	    double a[3], aa;
	    FOR_I3 a[i] = vert[(l+2)%3]->x[i] - vert[(l+1)%3]->x[i];
	    aa = m_ddot(3, a, a);

	    FOR_I3 {
	        D2tmp[i][i] += 0.5*aa;
		FOR_J3 D2tmp[i][j] -= 0.5*a[i]*a[j];
	    }
	} 
	else { 
	    int n = 3 - l - m;
	    double a1[3], a2[3], a1a2;
	    FOR_I3 {
	        a1[i] = vert[l]->x[i] - vert[n]->x[i];
		a2[i] = vert[m]->x[i] - vert[n]->x[i];
	    }
	    a1a2 = m_ddot(3, a1, a2);

	    FOR_I3 {
	        D2tmp[i][i] -= 0.5*a1a2;
		FOR_J3 D2tmp[i][j] += a1[i]*a2[j] - 0.5*a2[i]*a1[j];
	    }
	}

	// Now D2tmp = d^2 A^2/dx(l,:) dx(m,:).  We need to convert it
	// to d^2 A/dx[l] dx[m]
        FOR_I3
	FOR_J3 {
	    D2tmp[i][j] = -D1[l][i]*D1[m][j] + 0.5*D2tmp[i][j];
	}
	m_dscal(9, 1.0/area, *D2tmp);

	// Copy solution
	FOR_I3
	FOR_J3
	    D2[l][i][m][j] = D2tmp[i][j];

	if (l != m) {
	    FOR_I3
	    FOR_J3
		D2[m][j][l][i] = D2tmp[i][j];
	}
    }  // l, m
}


/* The minimum angle of a triangle */
double Tri::minAngle()
{
    // xx01 = x1 - x0
    // xx02 = x2 - x0
    // xx12 = x2 - x1
    double xx01[3], xx02[3], xx12[3];
    double len01, len02, len12;

    FOR_I3 {
        xx01[i] = vert[1]->x[i] - vert[0]->x[i];
        xx02[i] = vert[2]->x[i] - vert[0]->x[i];
        xx12[i] = vert[2]->x[i] - vert[1]->x[i];
    }
    
    len01 = normVec3D(xx01);
    len02 = normVec3D(xx02);
    len12 = normVec3D(xx12);

    double th0 = acos( m_ddot(3, xx01, xx02)/(len01*len02) );
    double th2 = acos( m_ddot(3, xx02, xx12)/(len02*len12) );
    double th1 = M_PI - th0 - th2;
    return min(th0, min(th1, th2));
}
