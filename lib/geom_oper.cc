#include "cxxheaders.h"
#include "geom_oper.h"
#include "mathfunc.h"
#include "marray.h"

/* double normVec3D
 * double normalizeVec3D
 * double distCart3D
 * void cross_product
 * void tri_normal
 * double tri_shapeFactor
 * double minDistToTri
 * bool point_in_seg
 * bool seg_overlap
 * double tri_solidAngle
 * void solidAngle
 * double bendAngle
 * void calcRotateMatrix
 * void rotate_vector
 * void calcPolygonAreaCenter
 * void slice_triangle 
 * bool tri_box_collide */

/* L2 norm a Cart3D vector */
double normVec3D(const double *x)
{
    return sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
}


/* Normalize a vector 
 * Arguments:
 *  x -- the vector to be normalized
 * Return the L2 norm of the original vector */
double normalizeVec3D(double *x)
{
    double r = sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
    double ir = 1.0/r;
    x[0] *= ir;
    x[1] *= ir;
    x[2] *= ir;
    return r;
}

// Distance between two points
double distCart3D(const double *a, const double *b)
{
    double x[3] = { a[0]-b[0], a[1]-b[1], a[2]-b[2] };
    return sqrt( x[0]*x[0] + x[1]*x[1] + x[2]*x[2] );
}

// cross product 
// c = a x b
void cross_product(const double *a, const double *b, double *c)
{
    c[0] = a[1]*b[2] - a[2]*b[1];
    c[1] = a[2]*b[0] - a[0]*b[2];
    c[2] = a[0]*b[1] - a[1]*b[0];
}


// triple product a * (b x c)
double triple_product(const double *a, const double *b, const double *c)
{
    return a[0]*(b[1]*c[2] - b[2]*c[1])
         + a[1]*(b[2]*c[0] - b[0]*c[2])
	 + a[2]*(b[0]*c[1] - b[1]*c[0]);
}


/* Calculate the normal and Jacobian of a triangle
 * Arguments:
 *  x0, x1, x2 -- vertex coordinates
 *  norm -- surface normal (unit vector)
 *  detj -- Jacobian, i.e. 2*area */
void tri_normal(const double *x0, const double *x1, const double *x2, 
		double *normal, double &detj)
{
    double a1[3], a2[3];

    for (int ii = 0; ii < 3; ++ii) {
	a1[ii] = x1[ii] - x0[ii];
	a2[ii] = x2[ii] - x0[ii];
    }

    cross_product(a1, a2, normal);
    detj = normalizeVec3D(normal);
}


void tri_normal(const double (*x)[3], double *normal, double &detj)
{
    tri_normal(x[0], x[1], x[2], normal, detj);
}


/* Triangle shape factor
 * Arguments:
 *  x0, x1, x2 -- vertex coordinates
 * Note:
 * 1. Defined as the longest edge/shortes edge
 */
double tri_shapeFactor(const double *x0, const double *x1, const double *x2)
{
    double L01 = 0.0, L12 = 0.0, L20 = 0.0;

    for (int ii = 0; ii < 3; ii++) {
        L01 += square(x1[ii] - x0[ii]);
        L12 += square(x2[ii] - x1[ii]);
	L20 += square(x0[ii] - x2[ii]);
    }

    double Lmax = sqrt( max(L01, max(L12, L20)) );
    double Lmin = sqrt( min(L01, min(L12, L20)) );
    return Lmax/Lmin;
}

/* The minimum distance from a point z to a triangle
 * Arguments:
 *  z -- the target point
 *  x[i] -- the i-th triangle corner point
 *  s, t -- the closest point is at x[0]*(1 - s - t) + x[1]*s + x[2]*t
 *
 * Reference:
 *  David Eberly, "Distance Between Point and Triangle in 3D" */
double minDistToTri(const double *z, const double (*x)[3], double &s, double &t)
{
    double dx1[3], dx2[3], dz[3], a, b, c, d, e, f, det, invDet;
    int whichRegion;

    for (int k = 0; k < 3; ++k) {
	dx1[k] = x[1][k] - x[0][k];
	dx2[k] = x[2][k] - x[0][k];
	dz[k] = x[0][k] - z[k];  // notice the order here
    }

    a = m_ddot(3, dx1, dx1);
    b = m_ddot(3, dx1, dx2);
    c = m_ddot(3, dx2, dx2);

    d = m_ddot(3, dx1, dz);
    e = m_ddot(3, dx2, dz);
    f = m_ddot(3, dz, dz);

    det = a*c - b*b; 
    invDet = 1.0/det;

    s = b*e - c*d;
    t = b*d - a*e; 

    if (s + t <= det) {
	if (s < 0)
	    whichRegion =  t < 0 ? 4 : 3;
	else 
	    whichRegion = t < 0 ? 5 : 0;
    } else {
	whichRegion = s < 0 ? 2 : (t < 0 ? 6 : 1);
    }

    switch (whichRegion) {
	case 2:
	    whichRegion = -(c+e) < 0 ? 3 : 1;
	    break;
	case 4:
	    whichRegion = d < 0 ? 5 : 3;
	    break;
	case 6:
	    whichRegion = -(a+d) < 0 ? 5 : 1;
	    break;
    }

    switch (whichRegion) {
      case 0:
	  s = invDet*s;
	  t = invDet*t;
	  break;

      case 1:
	  s = (c + e - b - d)/(a - 2*b + c);
	  s = min(1.0, max(0.0, s));
	  t = 1. - s;
	  break;

      case 3:
	  s = 0.;
	  t = -e/c;
	  t = min(1.0, max(0.0, t));
	  break;

      case 5:
	  t = 0.;
	  s = -d/a;
	  s = min(1.0, max(0.0, s));
	  break;
    }

    return sqrt(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f);
}


/* wheather a point lies in a line segment
 * Arguments:
 *   xmin, xmax -- end points of the segment
 *   x -- point coordinate
 *   L -- domain periodicity */
bool point_in_seg(double x, double xmin, double xmax, double L)
{
    if (x > xmin && x < xmax) return true;
    x -= floor((x - xmin)/L)*L;
    return (x < xmax);
}


/* Wheather two line segments intersect
 * Arguments:
 *   xmin, xmax -- end points of the 1st segment
 *   ymin, ymax -- end points of the 2nd segment
 *   L -- domain periodicity */
bool seg_overlap(double xmin, double xmax, double ymin, double ymax, double L)
{
    if (ymin >= xmin && ymin <= xmax) return true;
    if (ymax >= xmin && ymax <= xmax) return true;
    if (xmin >= ymin && xmin <= ymax) return true;
    if (xmax >= ymin && xmax <= ymax) return true;

    double dy = floor((ymin - xmin)/L)*L;
    ymin -= dy;
    ymax -= dy;
    return (ymin <= xmax || ymax >= xmin + L);
}


/* The solid angle formed between a point and a triangle
 * Arguments:
 *  x -- the target points
 *  y -- the triangle corner coordinates
 * Return: the solid angle
 *
 * Reference: Van Oosterom, A; Strackee, J (1983). "The Solid Angle of a Plane
 * Triangle". IEEE Trans. Biom. Eng. BME-30  (2): 125–126.  */
double tri_solidAngle(const double *x, const double (*y)[3])
{
    double a[3] = {y[0][0]-x[0], y[0][1]-x[1], y[0][2]-x[2]};
    double b[3] = {y[1][0]-x[0], y[1][1]-x[1], y[1][2]-x[2]};
    double c[3] = {y[2][0]-x[0], y[2][1]-x[1], y[2][2]-x[2]};

    double la = normVec3D(a);
    double lb = normVec3D(b);
    double lc = normVec3D(c);

    double abc = triple_product(a, b, c);
    double D = la*lb*lc + la*m_ddot(3,b,c) + lb*m_ddot(3,a,c) + lc*m_ddot(3,a,b);

    return 2*atan2(abc, D);
}


// Compute the inner solid angle at the corner of a simplex
// Arguments:
//   xc -- the corner point
//   N -- the number of neighboring points that surrounds xc
//   x -- coordinates of points surrounding xc, pointwise stored
//        i.e. x0,y0,z0, x1,y1,z1, ...
//   SA -- the solid angle
//   C -- the C-factor tensor, see notes
// Note:
//   -- The neighboring points surround xc in a clockwise fashion, i.e., the
//   normal of the triangule formed by xc, x(i), x(i+1) always points
//   outwards.
//
//   -- Reversing the order of x(i,:)s gives the supplementary solid angle, which
//   equals to (4 Pi - the original solid angle)
void solidAngle(const double *xc, int N, MArray<double,2> &x, double &SA, double (*C)[3])
{
    // axis(i,:) is the unit vector that points from xc to x(i,:)
    // normal(i,:) is the normal vector of the i-th face
    MArray<double,2> axis(N,3), normal(N,3);

    // initialize
    SA = 0;
    if (C) {
        for (int d0 = 0; d0 < 3; d0++)
	for (int d1 = 0; d1 < 3; d1++)
	    C[d0][d1] = 0;
    }

    // compute axis
    for (int i0 = 0; i0 < N; ++i0) {
        for (int d = 0; d < 3; ++d) {
	    axis(i0,d) = x(i0,d) - xc[d];
	}

	double s = m_dnrm2(3, &axis(i0,0));
	m_dscal(3, 1.0/s, &axis(i0,0));
    } // i0

    // compute normal of each side face
    // and the contribution to C factor
    for (int i0 = 0; i0 < N; ++i0) {
        int i1;
	double s, xctmp;

        i1 = modulo(i0+1, N);

	// surface normal
	cross_product(&axis(i0,0), &axis(i1,0), &normal(i0,0));
	s = m_dnrm2(3, &normal(i0,0));
	m_dscal(3, 1.0/s, &normal(i0,0));

	if (C) {
	    double s, th, xctmp[3];

	    s = m_ddot(3, &axis(i0,0), &axis(i1,0));
	    s = min(1.0, max(-1.0, s));
	    th = acos(s);

	    // xctmp = (3*PI/2)*(center of mass) * area 
	    // for the side wall
	    for (int d = 0; d < 3; ++d) {
		xctmp[d] = 0.5*(axis(i0,d) + axis(i1,d));
	    }
	    s = m_dnrm2(3, xctmp);
	    m_dscal(3, (1.0/M_PI)*sin(0.5*th)/s, xctmp);
            
	    for (int d0 = 0; d0 < 3; ++d0)
	    for (int d1 = 0; d1 < 3; ++d1) {
	        C[d0][d1] -= normal(i0,d0)*xctmp[d1];
	    } // d0, d1
	} // C
    } // i0


    // compute the solid angle
    for (int i0 = 0; i0 < N; ++i0) {
        int i1;
	double n01[3], costh, sinth;

	i1 = modulo(i0+1, N);

	// compute the rotational angle between normal(i0,:) and normal(i1,:)
	costh = m_ddot(3, &normal(i0,0), &normal(i1,0));

	cross_product(&normal(i0,0), &normal(i1,0), n01);
	sinth = -m_ddot(3, &axis(i1,0), n01);

	SA += M_PI - atan2(sinth, costh);
    }

    SA -= (N - 2)*M_PI;

    // add the isotropic part of the C-tensor
    if (C) {
        for (int d0 = 0; d0 < 3; ++d0) {
	    C[d0][d0] += SA/(2*M_PI);
	}
    }
}


/* Bending angle th
 * Arguments:
 *   x[0-3] -- the four points that form two triangles and 
 *             share a common edge
 *                   x1
 *                  /| \
 *                 / |  \
 *                /  |   \
 *               x2  |   x3
 *                \  |   /
 *                 \ |  /
 *                  \| /
 *                   x0
 *  -- The rotation of the two triangles are (012) and (031), they have
 *     normal vectors n0 and n1
 *  -- th is in [-pi, pi]
 *  -- costh = n0 * n1 
 *  -- sinth = (n1 - n0)*(x3 - x2) */
double bendAngle(const double *x0, const double *x1, const double *x2, const double *x3)
{
    double n0[3], n1[3], detj;
    double costh, dndx, th;

    tri_normal(x0, x1, x2, n0, detj);
    tri_normal(x0, x3, x1, n1, detj);
    costh = m_ddot(3, n0, n1);
    costh = min(1.0, max(-1.0, costh));

    dndx = 0.0;
    for (int ii = 0; ii < 3; ++ii)
        dndx += (n1[ii] - n0[ii])*(x3[ii] - x2[ii]);

    th = (dndx >= 0)? acos(costh) : -acos(costh);
    return th;
}


/* Calculate the rotation matrix that rotate xtri_old[i][:] to xtri_new[i][:] 
 * Arguments:
 *  xold[i] -- coordinates of the i-th vertex of the old triangle
 *  xnew[i] -- coordinates of the i-th vertex of the new triangle 
 *  A -- the rotation matrix
 * Note:
 *  -- Allow an additional rotation
 *  -- The old and new triangles must be of the same shape  */
void calcRotateMatrix(const double (*xold)[3], const double (*xnew)[3], double (*A)[3])
{
    // Solve A * B = C
    // where the columns of B are the two tangents and normal of the old triangle
    // the columns of C are those of the new triangle

    // We actually solve B^T A^T = C^T, i.e.
    //   A^T = B^{-T} C^T

    double B[3][3], C[3][3], iB[3][3];

    FOR_I3 {
        B[0][i] = xold[1][i] - xold[0][i];
	B[1][i] = xold[2][i] - xold[0][i];

	C[0][i] = xnew[1][i] - xnew[0][i];
	C[1][i] = xnew[2][i] - xnew[0][i];
    }

    cross_product(B[0], B[1], B[2]);
    normalizeVec3D(B[2]);

    cross_product(C[0], C[1], C[2]);
    normalizeVec3D(C[2]);

    invMat3(B, iB);	// iB = B^{-T}

    // A^T = B^{-T} * C^T
    FOR_I3
    FOR_J3 {
        A[j][i] = iB[i][0]*C[0][j] + iB[i][1]*C[1][j] + iB[i][2]*C[2][j];
    }
}


/* Calculate a rotation matrix
 * Arguments:
 *   oldx -- old x-axis
 *   newx -- new x-axis
 *   A -- rotating matrix, A*oldx = newx
 * Note:
 *   1. ||oldx|| = ||newx|| = 1.0 */
void calcRotateMatrix(const double *oldx, const double *newx, double (*A)[3])
{
    // Determine 2 vectors that are orthogonal to the z-axis 
    double z[3], oldy[3], newy[3];
    cross_product(oldx, newx, z);

    double rr = normalizeVec3D(z);

    if (rr < 1.E-10) {
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++) {
	    A[i][j] = (i == j)? 1.0 : 0.0;
	}
	return;
    } 
    else {
	cross_product(z, oldx, oldy);
	cross_product(z, newx, newy);

	FOR_I3
	FOR_J3 {
	    A[i][j] = newx[i]*oldx[j] + newy[i]*oldy[j] + z[i]*z[j];
	}
    }
}


/* Rotate a vector along some axis by a certain angle 
 * Arguments:
 *   axis -- rotate aixs
 *   th -- rotate angle
 *   x -- the vector to be rotated */
void rotate_vector(const double *axis, double th, double *x)
{
    // x0 parallel to axis, x1 perpendicular to axis
    double x0[3], x1[3], x2[3];	
    double costh=cos(th), sinth=sin(th);

    // Decompose the vector to components that are parallel and
    // perpendicular to the axis
    double c = m_ddot(3, axis, x);
    FOR_I3 {
        x0[i] = c*axis[i];
        x1[i] = x[i] - x0[i];
    }
    cross_product(axis, x1, x2);

    // Rotate the vector
    FOR_I3 x[i] = x0[i] + costh*x1[i] + sinth*x2[i];
}


/* Calculate the area and center of a polygon
 * Argument:
 *   np -- number of vertices
 *   xp[i] -- the i-th vertex
 *   area -- area
 *   center -- center
 * Note:
 *   1. Assume all vertices are co-planar
 *   2. Assume the polygon is convex */
void calcPolygonAreaCenter(int np, const double (*xp)[3], double &area, double *center)
{
    // Temporary center of the polygon
    double x0[3] = {0.0, 0.0, 0.0};
    for (int ip = 0; ip < np; ip++) {
	for (int ii = 0; ii < 3; ii++)
	  x0[ii] += xp[ip][ii];
    }
    for (int ii = 0; ii < 3; ii++) x0[ii] /= np;

    // Initialize
    area = 0.0;
    if (center) {
        for (int ii = 0; ii < 3; ii++) center[ii] = 0;
    }

    for (int i0 = 0; i0 < np; i0++) {
        int i1 = (i0+1)%np;
	double a0[3], a1[3];
	for (int ii = 0; ii < 3; ii++) {
	    a0[ii] = xp[i0][ii] - x0[ii];
	    a1[ii] = xp[i1][ii] - x0[ii];
	}

	double A = m_ddot(3, a0, a0);
	double B = m_ddot(3, a0, a1);
	double C = m_ddot(3, a1, a1);
	double areaTmp = 0.5*sqrt( fabs(A*C - B*B) );
	area += areaTmp;

	if (center) {
	    double xcTmp[3] = { 0.0, 0.0, 0.0 }; 
	    for (int ii = 0; ii < 3; ii++)
		center[ii] += areaTmp*THRD*(x0[ii] + xp[i0][ii] + xp[i1][ii]);
        }
    } // i0

    if (center) {
	if (area > 1.E-10) {
	    FOR_I3 center[i] /= area;
	} 
	else {
	    FOR_I3 center[i] = x0[i];
        }
    }
}

/* Slide a triangle using coordinate planse
 * Arguments:
 *   xtri[i] -- the i-th vertex coordinate
 *   normal -- normal direction of the planes
 *   zcut -- cutting plane coordinates
 *   nsub -- number of sub-polygons
 *   npsub -- number of vertices of sub-polygons
 *   xsub -- sub-polygon vertex coordinate
 * Note:
 *   The length of npsub >= zcut.size() + 1 
 */   
void slice_triangle(const double (*xtri)[3], const double normal[3], const vector<double> &zcut, 
	int &nsub, int *npsub, double (*xsub)[5][3])
{
    // Find the z-coordinate of each polygon vertex
    double ztri[3];
    int iztri[3];

    for (int i = 0; i < 3; i++) {
        ztri[i] = m_ddot(3, xtri[i], normal);

	vector<double>::const_iterator low;
	low = lower_bound(zcut.begin(), zcut.end(), ztri[i]);;
	iztri[i] = (int)(low - zcut.begin()) - 1;
    } // i

    int izMin = iztri[0], izMax = iztri[0];
    for (int i = 1; i < 3; i++) {
        izMin = min(izMin, iztri[i]);
	izMax = max(izMax, iztri[i]);
    }

    // Initialize
    nsub = izMax - izMin + 1;
    for (int iz = 0; iz < nsub; iz++) npsub[iz] = 0;

    for (int i0 = 0; i0 < 3; i0++) {
        int i1 = (i0 + 1)%3;
	int iz0 = iztri[i0];
	int iz1 = iztri[i1];

	int isub = iz0 - izMin;
	m_dcopy(3, xtri[i0], xsub[isub][npsub[isub]]);
	npsub[isub]++;

	// Scan-line
	int diz, iz_f, iz_l;
	if (iz1 == iz0) {
	    continue;
	} else if (iz1 > iz0) {
	    diz = 1;
	    iz_f = iz0 + 1;
	    iz_l = iz1 + 1;
	} else {
	    diz = -1;
	    iz_f = iz0;
	    iz_l = iz1;
	}


	for (int iz = iz_f; iz != iz_l; iz += diz) {
	    double s = (zcut[iz] - ztri[i0])/(ztri[i1] - ztri[i0]);
	    s = min(1.0, max(0.0, s));

	    // Create a new points at the intersection which belongs to two
	    // sub-polygons
	    double xnew[3];
	    for (int ii = 0; ii < 3; ii++) {
	        xnew[ii] = xtri[i0][ii] + s*(xtri[i1][ii] - xtri[i0][ii]);
	    }

	    isub = iz - izMin;
	    m_dcopy(3, xnew, xsub[isub][npsub[isub]]);
	    npsub[isub]++;

	    isub--;
	    m_dcopy(3, xnew, xsub[isub][npsub[isub]]);
	    npsub[isub]++;
	}
    }
}


/* Detect wheather a triangle collide with a box
 * Arguments:
 *   x[i] -- the i-th triangle vertex
 *   xmin, xmax -- the bounding box */
bool tri_box_collide(const double (*xtri)[3], double *xmin, double *xmax)
{
    // Translate the center of the box to origin
    double x[3][3], L[3];
    FOR_I3 {
        double c = 0.5*(xmin[i] + xmax[i]);
        x[0][i] = xtri[0][i] - c;
        x[1][i] = xtri[1][i] - c;
        x[2][i] = xtri[2][i] - c;

	L[i] = fabs(xmax[i] - xmin[i]);
    }

    // For the three axis
    double f[3][3], nrml[3];
    FOR_I3 {
        f[0][i] = x[0][i] - x[2][i];
        f[1][i] = x[1][i] - x[0][i];
        f[2][i] = x[2][i] - x[1][i];
    }
    cross_product( f[0], f[1], nrml );

    // Test against three axes
    FOR_I3 {
        double pmin = min( x[0][i], min(x[1][i], x[2][i]) );
	double pmax = max( x[0][i], max(x[1][i], x[2][i]) );
      
        if ( pmin > 0.5*L[i] || pmax < -0.5*L[i] ) return false;
    }


    // Test against triangle normal
    double r = 0.5*( fabs(nrml[0])*L[0] 
                   + fabs(nrml[1])*L[1] 
		   + fabs(nrml[2])*L[2] );
    double p = m_ddot(3, nrml, x[0]);
    if (p > r || p < -r) return false;

    // Test aginst the 9 vectors
    for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
        double a[3] = {0.0};
	a[i] = 1.0;

	double h[3];
	cross_product(a, f[j], h);

	double r = 0.5*( fabs(h[0])*L[0] 
                       + fabs(h[1])*L[1] 
		       + fabs(h[2])*L[2] );

	double p0 = m_ddot(3, h, x[0]);
	double p1 = m_ddot(3, h, x[1]);
	double p2 = m_ddot(3, h, x[2]);
	double pmin = min(p0, min(p1, p2));
	double pmax = max(p0, max(p1, p2));

	if (pmin > r || pmax < -r) return false;
    }

    return true;
}
