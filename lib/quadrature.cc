#include "cxxheaders.h"
#include "quadrature.h"
#include "mathfunc.h"
#include "tri.h"

/* quadrature::_init 
 * quadrature::select_rule_2d 
 * gauleg 
 * calcDuffyQuadPoint 
 * calcLogQuadPoint */

// Global variables
namespace quadrature
{
    int _inited = false; 
    Quad2D _Q_TRI_1;
    Quad2D _Q_TRI_3;
    Quad2D _Q_TRI_7;
    Quad2D _Q_QUAD_3;
    Quad2D _Q_QUAD_4;
}


// Initialize all quadratures used
void quadrature::_init()
{
    { // 1-point quadrature for triangle
	Quad2D &q = _Q_TRI_1;

	q.allocateMemory(1);

	q.w(0) = 0.5;
	q.x(0) = 1.0/3;
	q.y(0) = 1.0/3;
    }

    {    // 3-point quadrature for triangle
	double w, r, s, t;
	Quad2D &q = _Q_TRI_3;

	q.allocateMemory(3);

	r = 2./3;    s = t = 1./6; 
	w = 0.5/3;

	q.w(0) = q.w(1) = q.w(2) = w;
	q.x(0) = r;    q.y(0) = s;
	q.x(1) = s;    q.y(1) = t;
	q.x(2) = t;    q.y(2) = r;
    }


    {    // 7-point quadrature for triangle
	double w, r, s, t;
        Quad2D &q = _Q_TRI_7;

	q.allocateMemory(7);


	r = s = (6. - sqrt(15.))/21;    t = 1 - r - s;
	w = (155. - sqrt(15.))/2400;

	q.w(0) = q.w(1) = q.w(2) = w;

	q.x(0) = r;    q.y(0) = s;
	q.x(1) = s;    q.y(1) = t;
	q.x(2) = t;    q.y(2) = r;

	r = s = (6. + sqrt(15.))/21;    t = 1 - r - s;
	w = (155. + sqrt(15.))/2400;

	q.w(3) = q.w(4) = q.w(5) = w;

	q.x(3) = r;    q.y(3) = s;
	q.x(4) = s;    q.y(4) = t;
	q.x(5) = t;    q.y(5) = r;

	r = s = t = 1./3;
	w = 9./80;

	q.w(6) = w;
	q.x(6) = r;    q.y(6) = s;
    }


    {    // 3x3 point quadrilateral
	const int n = 3;
        double x[n], w[n];

	Quad2D &q = _Q_QUAD_3;
	q.allocateMemory(n*n);


	gauleg(-1.0, 1.0, n, x, w);
	int m = 0;
	for (int i = 0; i < n; ++i)
	for (int j = 0; j < n; ++j) {
	    q.w(m) = w[i]*w[j];
	    q.x(m) = x[i];
	    q.y(m) = x[j];
	    ++m;
	}
    }

    {    // 4x4 point quadrilateral
        const int n = 4;
	double x[n], w[n];

	Quad2D &q = _Q_QUAD_4;
	q.allocateMemory(n*n);

	gauleg(-1.0, 1.0, n, x, w);
	int m = 0;
	for (int i = 0; i < n; ++i)
	for (int j = 0; j < n; ++j) {
	    q.w(m) = w[i]*w[j];
	    q.x(m) = x[i];
	    q.y(m) = x[j];
	    ++m;
	}
    }
}


/* Select a quadrature rule
 * Arguments:
 *   ch -- name of the quadrature rule
 *         can be "TRI_3", "QUAD_3", "QUAD_4" */
Quad2D& quadrature::select_rule_2d(const char *ch)
{
    if (!_inited) {
        _init();
	_inited = true;
    }

    std::string qname(ch);

    if (qname.compare("TRI_1") == 0) {
        return _Q_TRI_1;
    } else if (qname.compare("TRI_3") == 0) {
        return _Q_TRI_3;
    } else if (qname.compare("TRI_7") == 0) {
        return _Q_TRI_7;
    } else if (qname.compare("QUAD_3") == 0) {
        return _Q_QUAD_3;
    } else if (qname.compare("QUAD_4") == 0) {
        return _Q_QUAD_4;
    } else {
        printf("Error: ");
        printf("can not find quadrature rule of name %s\n", ch);
	exit(0);
    }
}


/* Gauss-Legendre quadrature
 * Arguments:
 * [x1, x2] -- the interval
 * n -- the number of points
 * x, w -- the abscissas and weights */
int gauleg(double x1, double x2, int n, double *x, double *w)
{
    assert(n > 0);

    // Direct computation for small n
    if (n <= 8) {
	switch (n) {
	case 1:
	    x[0] = 0.0;				w[0] = 2.0;
	    break;
	case 2:
	    x[0] = 0.5773502691896257645091488;	w[0] = 1.0;
	    x[1] = -x[0];			w[1] = w[0];
	    break;
	case 3:
	    x[0] = 0;				w[0] = 0.8888888888888888888888889;
	    x[1] = 0.7745966692414833770358531;	w[1] = 0.5555555555555555555555556;
	    x[2] = -x[1];			w[2] = w[1];
	    break;
	case 4:
	    x[0] = 0.3399810435848562648026658;	w[0] = 0.6521451548625461426269361;
	    x[1] = -x[0];			w[1] = w[0];
	    x[2] = 0.8611363115940525752239465;	w[2] = 0.3478548451374538573730639;
	    x[3] = -x[2];			w[3] = w[2];
	    break;
	case 5:
	    x[0] = 0.0;				w[0] = 0.5688888888888888888888889;
	    x[1] = 0.5384693101056830910363144;	w[1] = 0.4786286704993664680412915;
	    x[2] = -x[1];			w[2] = w[1];
	    x[3] = 0.9061798459386639927976269;	w[3] = 0.2369268850561890875142640;
	    x[4] = -x[3];			w[4] = w[3];
	    break;
	case 6:
	    x[0] = 0.2386191860831969086305017;	w[0] = 0.4679139345726910473898703;
	    x[1] = -x[0];			w[1] = w[0];
	    x[2] = 0.6612093864662645136613996;	w[2] = 0.3607615730481386075698335;
	    x[3] = -x[2];			w[3] = w[2];
	    x[4] = 0.9324695142031520278123016;	w[4] = 0.1713244923791703450402961;
	    x[5] = -x[4];			w[5] = w[4];
	    break;
	case 7:
	    x[0] = 0.0;				w[0] = 0.4179591836734693877551020;
	    x[1] = 0.4058451513773971669066064;	w[1] = 0.3818300505051189449503698;
	    x[2] = -x[1];			w[2] = w[1];
	    x[3] = 0.7415311855993944398638648;	w[3] = 0.2797053914892766679014678;
	    x[4] = -x[3];			w[4] = w[3];
	    x[5] = 0.9491079123427585245261897;	w[5] = 0.1294849661688696932706114;
	    x[6] = -x[5];			w[6] = w[5];
	    break;
	case 8:
	    x[0] = 0.1834346424956498049394761;	w[0] = 0.3626837833783619829651504;
	    x[1] = -x[0];			w[1] = w[0];
	    x[2] = 0.5255324099163289858177390;	w[2] = 0.3137066458778872873379622;
	    x[3] = -x[2];			w[3] = w[2];
	    x[4] = 0.7966664774136267395915539;	w[4] = 0.2223810344533744705443560;
	    x[5] = -x[4];			w[5] = w[4];
	    x[6] = 0.9602898564975362316835609;	w[6] = 0.1012285362903762591525314;
	    x[7] = -x[6];			w[7] = w[6];
	    break;
	}

	double xm = 0.5*(x2 + x1);
	double xl = 0.5*(x2 - x1);
	for (int i = 0; i < n; i++) x[i] = xm + xl*x[i];
	for (int i = 0; i < n; i++) w[i] *= xl;

	return 0;
    }

    const double PI = M_PI;
    const double EPS=1.0e-14;
    double z1, z, xm, xl, pp, p3, p2, p1;
    int m = (n+1)/2;

    xm = 0.5*(x2 + x1);
    xl = 0.5*(x2 - x1);

    for (int i=0;i<m;i++) {
	z=cos(PI*(i+0.75)/(n+0.5));
	do {
	    p1=1.0;
	    p2=0.0;
	    for (int j=0;j<n;j++) {
		p3=p2;
		p2=p1;
		p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
	    }
	    pp=n*(z*p1-p2)/(z*z-1.0);
	    z1=z;
	    z=z1-p1/pp;
	} while (fabs(z-z1) > EPS);

	x[i]=xm-xl*z;
	x[n-1-i]=xm+xl*z;
	w[i]=2.0*xl/((1.0-z*z)*pp*pp);
	w[n-1-i]=w[i];
    }

    return 0;
}


/* Calculate Duffy quadrature rule
 * Arguments:
 *   n -- number of quadrature points in each direction
 *   sq -- s-coord of quadrature points
 *   tq -- t-coord of quadrature points 
 *   wq -- Jacobian (for a unit triangle of area 1/2) 
 * Note:
 *   -- The reference triangle is (0,0), (0,1), and (1,0) */
void calcDuffyQuadPoint(int n, vector<double> &sq, vector<double> &tq, vector<double> &wq)
{
    assert(n > 0);

    static int n_save = -1;
    static vector<double> sq_save, tq_save, wq_save;

    // Re-use saved quadrature points when applicable
    if (n_save == n) {
        sq = sq_save;
	tq = tq_save;
	wq = wq_save;
	return;
    }

    sq.resize(n*n);
    tq.resize(n*n);
    wq.resize(n*n);

    double *xq1d = new double[n];
    double *wq1d = new double[n];
    gauleg(0, 1, n, xq1d, wq1d);

    for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
        int p = i*n + j;
        sq[p] = xq1d[i];
	tq[p] = xq1d[i]*xq1d[j];
	wq[p] = xq1d[i]*wq1d[i]*wq1d[j];

	// Transform the third point from (1,1) to (0,1)
	sq[p] -= tq[p];		
    }

    delete [] xq1d;
    delete [] wq1d;

    // Save result for re-use
    n_save = n;
    sq_save = sq;
    tq_save = tq;
    wq_save = wq;
}


/* Calculate Duffy quadrature rule
 * Arguments:
 *   s0, t0 -- coord of quadrature points 
 *   n -- number of quadrature points in each direction
 *   sq -- s-coord of quadrature points
 *   tq -- t-coord of quadrature points 
 *   wq -- quadrature weight */
void calcDuffyQuadPoint(double s0, double t0, int n, 
		vector<double> &sq, vector<double> &tq, vector<double> &wq)
{
    sq.reserve(3*n*n);
    tq.reserve(3*n*n);
    wq.reserve(3*n*n);

    sq.resize(0);
    tq.resize(0);
    wq.resize(0);

    for (int isub = 0; isub < 3; ++isub) {
        double detjSub, s1, s2, t1, t2;
        // (s0,t0), (s1,t1), (s2,t2) are the reference coordinates of
        //  the three corners of the sub-triangle in the original triangle
        switch (isub) {
        case 0:
            s1 = 1.0;    t1 = 0.0;
            s2 = 0.0;    t2 = 1.0;
            detjSub = (1 - s0 - t0);
            break;
        case 1:
            s1 = 0.0;    t1 = 1.0;
            s2 = 0.0;    t2 = 0.0;
            detjSub = s0;
            break;
        case 2:
            s1 = 0.0;    t1 = 0.0;
            s2 = 1.0;    t2 = 0.0;
            detjSub = t0;
            break;
        }
        if (fabs(detjSub) < 1.E-10) continue;

	vector<double> sq_tmp, tq_tmp, wq_tmp;
	calcDuffyQuadPoint(n, sq_tmp, tq_tmp, wq_tmp);
	for (int i = 0; i < sq_tmp.size(); i++) {
	    double sloc = sq_tmp[i];
	    double tloc = tq_tmp[i];

            // Map: (0,0) -> (s0,t0)
	    //      (1,0) -> (s1,t1)
	    //      (0,1) -> (s2,t2)
	    sq_tmp[i] = s0 + (s1-s0)*sloc + (s2-s0)*tloc;
	    tq_tmp[i] = t0 + (t1-t0)*sloc + (t2-t0)*tloc;
	    wq_tmp[i] *= detjSub;
	}

	sq.insert(sq.end(), sq_tmp.begin(), sq_tmp.end());
	tq.insert(tq.end(), tq_tmp.begin(), tq_tmp.end());
	wq.insert(wq.end(), wq_tmp.begin(), wq_tmp.end());
    }
}

/* Calcualte the quadrature points using Log transform to remove
 * nearly singularities
 * Arguments:
 *   dist -- the distance from the target to the plane
 *   h -- the distance from the project point to the edge
 *   thmin, thmax -- the min and max angles
 *   n -- the number of quadrature points in each direction
 *   sq, tq, dA -- quadrature coordinates and area
 */
void calcLogQuadPoint(double dist, double h, double thmin, double thmax,
		     int n, vector<double> &sq, vector<double> &tq, vector<double> &dA)
{
   assert(n > 0);

   sq.resize(n*n);
   tq.resize(n*n);
   dA.resize(n*n);

   double A[2][2];
   A[0][0] = h;	A[1][0] = h*tan(thmin);
   A[0][1] = h;	A[1][1] = h*tan(thmax);

   double invA[2][2];
   invMat2(A, invA);

   double *t = new double[n];
   double *wt = new double[n];
   double *R = new double[n];
   double *wR = new double[n];

   double tmin = 0.5*h*log( (1+sin(thmin))/(1 - sin(thmin)) );
   double tmax = 0.5*h*log( (1+sin(thmax))/(1 - sin(thmax)) );
   gauleg(tmin, tmax, n, t, wt);

   for (int i = 0; i < n; i++) {
       double sinth = (exp(2*t[i]/h) - 1)/(exp(2*t[i]/h) + 1);
       double costh = sqrt(1.0 - sinth*sinth);

       double rhomin = 0.0;
       double rhomax = h/costh;

       double Rmin = log(rhomin + dist);
       double Rmax = log(rhomax + dist);
       gauleg(Rmin, Rmax, n, R, wR);

       for (int j = 0; j < n; j++) {
           double rho = exp(R[j]) - dist;
           double x = rho*costh;
           double y = rho*sinth;

           int p = i*n + j;
           sq[p] = invA[0][0]*x + invA[0][1]*y;
           tq[p] = invA[1][0]*x + invA[1][1]*y;
           dA[p] = wt[i]*wR[j]*rho*costh/h*(rho + dist);
       }
   }

   delete [] t;
   delete [] wt;
   delete [] R;
   delete [] wR;
}


/* Calcualte the quadrature points using Log transform to remove
 * nearly singularities
 * Arguments:
 *   xtar -- the observer point
 *   T -- the triangle
 *   n -- the number of quadrature points in each direction
 *   sq, tq, wq -- quadrature coordinates and weights
 */
void calcLogQuadPoint(const double *xtar, const Tri &T,
		     int n, vector<double> &sq, vector<double> &tq, vector<double> &dA)
{
    sq.reserve(3*n*n);
    tq.reserve(3*n*n);
    dA.reserve(3*n*n);

    sq.resize(0);
    tq.resize(0);
    dA.resize(0);

    // Find the projection in the triangle
    double s0, t0;
    T.baryCoord(xtar, s0, t0);

    double x0[3];
    FOR_I3 x0[i] = (1 - s0 - t0)*T.vert[0]->x[i]
                + s0*T.vert[1]->x[i]
		+ t0*T.vert[2]->x[i];

    // The distance to the triangle
    double dist = 0.0;
    FOR_I3 dist += (xtar[i] - x0[i])*T.normal[i];

    for (int isub = 0; isub < 3; isub++) {
        double detjSub, s1, s2, t1, t2;
        // (s0,t0), (s1,t1), (s2,t2) are the reference coordinates of
        //  the three corners of the sub-triangle in the original triangle
        switch (isub) {
        case 0:
            s1 = 1.0;    t1 = 0.0;
            s2 = 0.0;    t2 = 1.0;
            detjSub = (1 - s0 - t0);
            break;
        case 1:
            s1 = 0.0;    t1 = 1.0;
            s2 = 0.0;    t2 = 0.0;
            detjSub = s0;
            break;
        case 2:
            s1 = 0.0;    t1 = 0.0;
            s2 = 1.0;    t2 = 0.0;
            detjSub = t0;
            break;
        }
        if (fabs(detjSub) < 1.E-6) continue;

        // Find the distance from x0 to the edge
        int i1 = (isub + 1)%3;
        int i2 = (isub + 2)%3;

        const double *x1 = T.vert[i1]->x;
        const double *x2 = T.vert[i2]->x;

        double x01[3];
        FOR_D3 x01[d] = x0[d] - x1[d];

        double x21[3], len21;
        FOR_D3 x21[d] = x2[d] - x1[d];
        len21 = sqrt(m_ddot(3, x21, x21));

        // scut is the projection of x0 on the edge
        double scut = m_ddot(3, x01, x21)/m_ddot(3, x21, x21);
        double xcut[3];
        FOR_D3 xcut[d] = x1[d] + scut*(x2[d] - x1[d]);

        double xx[3];
        FOR_D3 xx[d] = x0[d] - xcut[d];
        double h = sqrt(m_ddot(3, xx, xx));

        double thmin = atan(-scut*len21/h);
        double thmax = atan((1 - scut)*len21/h);

	vector<double> sq_tmp, tq_tmp, dA_tmp;
	calcLogQuadPoint(fabs(dist), h, thmin, thmax, n, sq_tmp, tq_tmp, dA_tmp);
	for (int i = 0; i < sq_tmp.size(); i++) {
	    double sloc = sq_tmp[i];
	    double tloc = tq_tmp[i];

            // Map: (0,0) -> (s0,t0)
	    //      (1,0) -> (s1,t1)
	    //      (0,1) -> (s2,t2)
	    sq_tmp[i] = s0 + (s1-s0)*sloc + (s2-s0)*tloc;
	    tq_tmp[i] = t0 + (t1-t0)*sloc + (t2-t0)*tloc;

	    if (detjSub < 0) dA_tmp[i] = -dA_tmp[i];
	}

	sq.insert(sq.end(), sq_tmp.begin(), sq_tmp.end());
	tq.insert(tq.end(), tq_tmp.begin(), tq_tmp.end());
	dA.insert(dA.end(), dA_tmp.begin(), dA_tmp.end());
    }
}
