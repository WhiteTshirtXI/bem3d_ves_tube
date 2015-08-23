#include "cxxheaders.h"
#include "vesicle.h"
#include "edge.h"
#include "membrane.h"
#include "mathfunc.h"
#include "mypetsc.h"
#include "mpi.h"

/* Vesicle::Vesicle 
 * Vesicle::~Vesicle 
 * Vesicle::buildControlPoints
 * Vesicle::calcQuadPoint
 * Vesicle::bendEnergy
 * Vesicle::bendForce
 * Vesicle::diffBendForce 
 * Vesicle::calcMeshRelaxVel
 * Vesicle::relaxMesh
 * Vesicle::gradVolArea
 * Vesicle::adjustVolArea 
 * Vesicle::vesaddtension (new function)*/

Vesicle::Vesicle() : Mesh()
{ 
}


Vesicle::~Vesicle()
{
}

/* Rebuild vertex valences and build control point list
 * Note:
 *  -- This needs to be done once and only once when the mesh
 *     topology is changed */
void Vesicle::buildControlPoints()
{
    int nvert = numVerts();
    int nface = numFaces();
    int (*f2v)[3] = new int[nface][3];
    getConnectivities(f2v);

    // Compute vertex valences
    vertValence.resize(nvert);
    std::fill(vertValence.begin(), vertValence.end(), 0);

    for (int iface = 0; iface < nface; iface++) {
        for (int l = 0; l < 3; l++) {
            int ivert = f2v[iface][l];
            vertValence[ivert]++;
        }
    }

    // Build control points
    subdiv::buildControlPoints(nvert, nface, f2v, ctrlPts);

    delete [] f2v;
}

void Vesicle::calcQuadPoint(
	const DoubArray1d &f, const DoubArray2d &df, const DoubArray3d &ddf, 
	const DoubArray2d &xctrl, QuadPoint &q)
{
    const int N = f.size(0);
    // x, dx, ddx, normal
    m_dclear(3, q.x);
    m_dclear(6, *q.dx);
    m_dclear(12, **q.ddx);

    for (int i = 0; i < N; i++) {
        const double *xi = &xctrl(i,0);
	m_daxpy(3, f(i), xi, q.x);

	m_daxpy(3, df(i,0), xi, q.dx[0]);
	m_daxpy(3, df(i,1), xi, q.dx[1]);

	m_daxpy(3, ddf(i,0,0), xi, q.ddx[0][0]);
	m_daxpy(3, ddf(i,1,1), xi, q.ddx[1][1]);
	m_daxpy(3, ddf(i,0,1), xi, q.ddx[0][1]);
    }
    m_dcopy(3, q.ddx[0][1], q.ddx[1][0]);

    cross_product(q.dx[0], q.dx[1], q.nrml);
    normalizeVec3D(q.nrml);

    // Metric tensor
    q.a[0][0] = m_ddot(3, q.dx[0], q.dx[0]);
    q.a[1][1] = m_ddot(3, q.dx[1], q.dx[1]);
    q.a[0][1] = q.a[1][0] = m_ddot(3, q.dx[0], q.dx[1]);

    // Reciprocal vectors and tensor
    invMat2(q.a, q.a_rcp);

    FOR_J3 {
        q.dx_rcp[0][j] = q.a_rcp[0][0]*q.dx[0][j] + q.a_rcp[0][1]*q.dx[1][j];
        q.dx_rcp[1][j] = q.a_rcp[1][0]*q.dx[0][j] + q.a_rcp[1][1]*q.dx[1][j];
    }

    // Curvature tensor
    q.b[0][0] = m_ddot(3, q.ddx[0][0], q.nrml);
    q.b[1][1] = m_ddot(3, q.ddx[1][1], q.nrml);
    q.b[0][1] = q.b[1][0] = m_ddot(3, q.ddx[0][1], q.nrml);

    for (int i = 0; i < 2; i++)
    for (int j = 0; j < 2; j++) {
	q.b_rcp[i][j] = q.a_rcp[i][0]*q.a_rcp[j][0]*q.b[0][0]
                      + q.a_rcp[i][1]*q.a_rcp[j][1]*q.b[1][1]
	              + q.a_rcp[i][0]*q.a_rcp[j][1]*q.b[0][1]
		      + q.a_rcp[i][1]*q.a_rcp[j][0]*q.b[1][0];
    }

    double detA = detMat2(q.a);
    q.area = sqrt(detA)*q.wght;
    q.H = 0.5*m_ddot(4, *q.a_rcp, *q.b);
    q.K = detMat2(q.b)/detA;
}


double Vesicle::bendEnergy()
{
    // Init
    double E = 0.0;

    // Quadrature rule
    const Quad2D &qrule = quadrature::select_rule_2d("TRI_3");
    QuadPoint q;

    for (int iface = 0; iface < numFaces(); iface++) {
        const Tri &face = faces[iface];

	int i0 = face.ivert[0];
	int i1 = face.ivert[1];
	int i2 = face.ivert[2];

	int N0 = vertValence[i0];
	int N1 = vertValence[i1];
	int N2 = vertValence[i2];

	int N = ctrlPts[iface].size();
	MArray<double,1> f(N);
	MArray<double,2> df(N,2);
	MArray<double,3> ddf(N,2,2);
	MArray<double,2> xctrl(N,3);

	MArray<double,1> lhs_loc(N), rhs_loc(N);
	lhs_loc = 0.0;
	rhs_loc = 0.0;

	for (int i = 0; i < N; i++) {
	    int ivert = ctrlPts[iface][i];
	    m_dcopy(3, verts[ivert].x, &xctrl(i,0));
	}


	for (int iq = 0; iq < qrule.n(); iq++) {
	    q.s = qrule.x(iq);
	    q.t = qrule.y(iq);
	    q.wght = qrule.w(iq);

	    subdiv::calcFunc(N0, N1, N2, q.s, q.t, f, df, ddf);
	    calcQuadPoint(f, df, ddf, xctrl, q);

	    E += 2*EB*q.H*q.H*q.area;
	} // iq
    } // iface

    return E;
}


/* Bending force 
 * Arguments:
 *  fb -- the bending force density at each vertex 
 * Note:
 *  Prerequisites: vertArea, qpoints */
void Vesicle::bendForce(MArray<double,2> &fb)
{
    int nvert = numVerts();
    int nface = numFaces();

    MArray<double,1> lhs(nvert);
    lhs = 0.0;

    // Init
    fb.resize(nvert,3);
    fb = 0.0;

    // First calculate f_a = dW/dx_a
    const Quad2D &qrule = quadrature::select_rule_2d("TRI_3");
    QuadPoint q;

    for (int iface = 0; iface < nface; iface++) {
        const Tri &face = faces[iface];

	int i0 = face.ivert[0];
	int i1 = face.ivert[1];
	int i2 = face.ivert[2];

	int N0 = vertValence[i0];
	int N1 = vertValence[i1];
	int N2 = vertValence[i2];

	int N = ctrlPts[iface].size();
	MArray<double,1> f(N);
	MArray<double,2> df(N,2);
	MArray<double,3> ddf(N,2,2);
	MArray<double,2> xctrl(N,3);

	MArray<double,1> lhs_loc(N);
	lhs_loc = 0.0;

	MArray<double,2> fb_loc(N,3);
	fb_loc = 0.0;

	for (int i = 0; i < N; i++) {
	    int ivert = ctrlPts[iface][i];
	    m_dcopy(3, verts[ivert].x, &xctrl(i,0));
	}


	for (int iq = 0; iq < qrule.n(); iq++) {
	    q.s = qrule.x(iq);
	    q.t = qrule.y(iq);
	    q.wght = qrule.w(iq);

	    subdiv::calcFunc(N0, N1, N2, q.s, q.t, f, df, ddf);
	    calcQuadPoint(f, df, ddf, xctrl, q);

	    // dH = 0.5*( A^i dot d_dx_i + ddx^{ij} n dot d_ddx_{ij} )
	    double A[2][3] = {0.0};
	    for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
		    m_daxpy(3, -2*q.b_rcp[i][j], q.dx[j], A[i]);
		}

		double c = 0.0;
		for (int j = 0; j < 2; j++)
		for (int k = 0; k < 2; k++) {
		    c += q.a_rcp[j][k]*m_ddot(3, q.dx_rcp[i], q.ddx[j][k]);
		}
		m_daxpy(3, -c, q.nrml, A[i]);
	    }	

	    for (int alph = 0; alph < N; alph++) {
		// d_H = dH/dx_{alpha}
	        double d_H[3] = {0.0};

		for (int i = 0; i < 2; i++) {
		    m_daxpy(3, 0.5*df(alph,i), A[i], d_H);
	        }

		double c = 0.0;
		for (int i = 0; i < 2; i++)
		for (int j = 0; j < 2; j++) {
		    c += q.a_rcp[i][j]*ddf(alph,i,j);
		}
		m_daxpy(3, 0.5*c, q.nrml, d_H);

		// d_J = 1/J dJ/dx_{alpha}, where J = sqrt(detA)
	        double d_J[3] = {0.0};
		for (int i = 0; i < 2; i++) {
		    m_daxpy(3, df(alph,i), q.dx_rcp[i], d_J);
		}

		m_daxpy(3, 2*EB*2*q.H*q.area, d_H, &fb_loc(alph,0));
		m_daxpy(3, 2*EB*q.H*q.H*q.area, d_J, &fb_loc(alph,0));

		lhs_loc(alph) += f(alph)*q.area;
	    } // alph
	} // iq

	// Add to global lhs and rhs
	for (int alph = 0; alph < N; alph++) {
	    int i = ctrlPts[iface][alph];

	    lhs(i) += lhs_loc(alph);
	    m_dadd(3, &fb_loc(alph,0), &fb(i,0));
	}
    } // iface

    for (int ivert = 0; ivert < nvert; ivert++) {
//	m_dscal(3, 1.0/lhs(ivert), &fb(ivert,0));
	m_dscal(3, 1.0/vertArea(ivert), &fb(ivert,0));
    }
} 


/* The change to bending force due to surface perturbation (linearized) 
 * Arguments:
 *  dx -- surface perturbation
 *  df -- change to bending force */
void Vesicle::diffBendForce(const MArray<double,2> &dx, MArray<double,2> &df)
{
    const int nvert = numVerts();
    const int nface = numFaces();

    // Save old coordinates
    MArray<double,2> xsave(nvert,3);
    getCoords(xsave);

    // Determine the scaling 
    // Average mesh size
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

    const double scal = 1.E-3*meshSize/max_dx;

    // Bend force after shape perturbation
    MArray<double,2> fbend_new(nvert,3);

    for (int ivert = 0; ivert < nvert; ivert++) {
        Point &vert = verts[ivert];
	m_daxpy(3, scal, &dx(ivert,0), vert.x);
    }
    bendForce(fbend_new);

    // Restore shape and bend force
    setCoords(xsave);
    updateGeometry();

    // Finite-difference
    df.resize(nvert,3);
    df = fbend_new;
    df -= fbend;
    df *= 1.0/scal;
}


/* Calculate the surface velocity for mesh relaxation
 * Arguments:
 *   v -- the relaxation velocity */
void Vesicle::calcMeshRelaxVel(MArray<double,2> &v)
{
    int nvert = numVerts();
    int nface = numFaces();
    int nedge = numEdges();

    // Init
    v.resize(nvert,3);
    v = 0.0;

    // Vertex valence
    MArray<int,1> valence(nvert);
    valence = 0;
    for (int iface = 0; iface < nface; iface++) {
        for (int l = 0; l < 3; l++) {
            int ivert = faces[iface].ivert[l];
            valence(ivert)++;
        }
    }

    // Calculate sprint constant
    MArray<double,1> sk(nvert);	

    sk = 1.0;	// Uniform sprint constant 

    // Weighted by the local curvature and 1/(local average triangle area)
    for (int ivert = 0; ivert < nvert; ivert++) {
	sk(ivert) *= pow(fabs(vertH(ivert)), 0.5);
	sk(ivert) *= pow(vertArea(ivert)/valence(ivert), 1.0);
    }

    // Normalize the sprint constant
    double mean_sk = m_dsum(nvert, sk.data())/nvert;
    sk *= 1./mean_sk;
    for (int ivert = 0; ivert < nvert; ivert++) {
        sk(ivert) = min(5.0, max(0.2, sk(ivert)));
    }

    // Calculate v
    MArray<double,1> dnrm(nvert);

    v = 0.0;
    dnrm = 0.0;

    for (int ied = 0; ied < nedge; ied++) {
	Edge &edge = edges[ied];

	int iv0 = edge.ivert[0];
	int iv1 = edge.ivert[1];

	double xx[3];
	for (int d = 0; d < 3; d++) xx[d] = verts[iv1].x[d] - verts[iv0].x[d];

	m_daxpy(3, sk(iv1), xx, &v(iv0,0));
	dnrm(iv0) += sk(iv1);

	m_daxpy(3, -sk(iv0), xx, &v(iv1,0));
	dnrm(iv1) += sk(iv0);
    } //ied

    for (int ivert = 0; ivert < nvert; ivert++) {
        m_dscal(3, 1./dnrm(ivert), &v(ivert,0)); 
    }

    // Project the moving vector onto the tangent plane
    for (int ivert = 0; ivert < nvert; ivert++) {
	Point &vert = verts[ivert];

	double vn = m_ddot(3, &v(ivert,0), &vertNrml(ivert,0));
	m_daxpy(3, -vn, &vertNrml(ivert,0), &v(ivert,0));
    } 
}


/* Relax the mesh
 * Return 0 if the mesh points are not moved
 *        1 otherwise 
 * Argument:
 *   niter -- number of relaxation iterations 
 *   str -- relaxation strength */
int Vesicle::relaxMesh(int niter, double str)
{
    const double MIN_ANGLE_THRESH = M_PI/6.0;
    int mesh_relaxed = 0;

    const int nvert = numVerts();
    const int nface = numFaces();

    // Relax the mesh
    for (int iter = 0; iter < niter; iter++) {
	updateGeometry();
	if (minAngle() > MIN_ANGLE_THRESH) break;

	MArray<double,2> v(nvert,3);
	calcMeshRelaxVel(v);

	for (int ivert = 0; ivert < nvert; ivert++) {
	    m_daxpy(3, str, &v(ivert,0), verts[ivert].x);
	}

	mesh_relaxed = 1;
    }

    return mesh_relaxed;
}


/* Compute the gradient of volume and area with respect to the
 * movement in normal direction
 * Arguments:
 *  gradV -- gradient of volume
 *  gradA -- gradient of area */
void Vesicle::gradVolArea(DoubArray1d &gradV, DoubArray1d &gradA)
{
    // Init
    int nvert = numVerts();
    gradV.resize(nvert);
    gradA.resize(nvert);

    gradV = 0.0;
    gradA = 0.0;

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &T = faces[iface];
	double dAdx[3][3];
	T.gradArea(dAdx);

        for (int l = 0; l < 3; l++) {
	    int ivert = T.ivert[l];
	    gradV(ivert) += THRD*m_ddot(3, &vertNrml(ivert,0), T.normal)*T.area;
	    gradA(ivert) += m_ddot(3, &vertNrml(ivert,0), dAdx[l]);
	}
    }
}


/* Adjust volume and area toward volTar and areaTar 
 * Arguments:
 *  Vtar -- target volume
 *  Atar -- target area */
void Vesicle::adjustVolArea(double Vtar, double Atar)
{
    int nvert = numVerts();
    MArray<double,1> gradV(nvert), gradA(nvert);
    gradVolArea(gradV, gradA);

    double lhs[2][2], ilhs[2][2];
    double rhs[2], sol[2];

    lhs[0][0] = m_ddot(nvert, gradV.data(), gradV.data());
    lhs[0][1] = m_ddot(nvert, gradV.data(), gradA.data());
    lhs[1][0] = lhs[0][1];
    lhs[1][1] = m_ddot(nvert, gradA.data(), gradA.data());

    rhs[0] = Vtar - vol;
    rhs[1] = Atar - area;

    invMat2(lhs, ilhs);

    sol[0] = ilhs[0][0]*rhs[0] + ilhs[0][1]*rhs[1];
    sol[1] = ilhs[1][0]*rhs[0] + ilhs[1][1]*rhs[1];

    for (int ivert = 0; ivert < nvert; ivert++) {
        double dn = sol[0]*gradV(ivert) + sol[1]*gradA(ivert);
	m_daxpy(3, dn, &vertNrml(ivert,0), verts[ivert].x);
    }
}

/* New function by spann 
 * for easier compatibility with stokesys calcstresslet functions
 *  analogous to cell.cc's calcSurfforce, but called after sigma is corrected and fbend is known
 *  Sets f=fbend+ftens*/
void Vesicle::vesaddtension()
{
	int nvert = numVerts();

    fbend.resize(nvert,3);
    ftens.resize(nvert,3);
    f.resize(nvert,3);

	//bendingForce already computed by now
	tensionForce(sigma, ftens);

    f = 0.0;
    f += fbend;
    f += ftens;
	
}
