#include "cxxheaders.h"
#include "redcell.h"
#include "mathfunc.h"
#include "membrane.h"
#include "subdiv.h"

/* RedCell::RedCell 
 * RedCell::~RedCell 
 * RedCell::updateGeometry
 * RedCell::tri_elasticity
 * RedCell::tri_elasticity_Skalak
 * RedCell::tri_elasticity_PNAS2002
 * RedCell::elasticForce 
 * RedCell::elasticEnergy 
 * RedCell::bendEnergy
 * RedCell::bendForce 
 * RedCell::tensionForce 
 * RedCell::viscousForce
 * RedCell::tweezerForce
 * RedCell::diffElasBendForce
 * RedCell::buildControlPoints
 * RedCell::calcQuadPoint */

RedCell::RedCell()
{ 
    cellRef = NULL;

    ES = 0.0;
    EB = 0.0;
    muM = 0.0;

    sigma = 0.0;
    area_E0 = 0.0;
    area_E1 = 0.0;
    area_E2 = 0.0;
}


RedCell::~RedCell()
{ 
}


void RedCell::updateGeometry()
{
    Mesh::updateGeometry();

    if (numEdges() == 0) buildEdgeList();	// needed by calcDoubleLayerJump
    calcDoubleLayerJump();
}


double RedCell::tri_elasticity(Tri &T, Tri &Tref, double (*grad)[3])
{
    switch (elas_model) { 
        case SKALAK_MODEL:
            return tri_elasticity_Skalak(T, Tref, grad);
	    break;
        case PNAS2002_MODEL:
            return tri_elasticity_PNAS2002(T, Tref, grad);
	    break;
        default:
	    throw(-1);
    }
}


double RedCell::tri_elasticity_Skalak(Tri &T, Tri &Tref, double (*grad)[3])
{
    // tangents
    double a0[3], a1[3];
    for (int j = 0; j < 3; ++j) {
	a0[j] = T.vert[0]->x[j] - T.vert[2]->x[j];
	a1[j] = T.vert[1]->x[j] - T.vert[2]->x[j];
    }

    // Strain invariants
    double I1 = m_ddot(4, &T.a[0][0], &Tref.arcp[0][0]) - 2; 
    double I2 = T.detA/Tref.detA - 1;

    if (grad) {
	m_dclear(9, *grad);

	double W1 = 0.5*ES*(1 + 0.5*I1);
	double W2 = -0.5*ES + 0.25*ED*I2;

	for (int i = 0; i < 2; ++i) {
	    double c = W1*2*Tref.arcp[i][0] + W2*2*(T.detA/Tref.detA)*T.arcp[i][0];
	    m_daxpy(3, c, a0, grad[i]);

	    c = W1*2*Tref.arcp[i][1] + W2*2*(T.detA/Tref.detA)*T.arcp[i][1];
	    m_daxpy(3, c, a1, grad[i]);
	} // i

	for (int j = 0; j < 3; ++j) {
	    grad[2][j] = -grad[0][j] - grad[1][j];
	}

	// Scale by triangle area
	m_dscal(9, Tref.area, *grad);
    }

    return Tref.area*( 0.5*ES*(I1 - I2 + 0.25*I1*I1) + 0.125*ED*I2*I2 );
}



double RedCell::tri_elasticity_PNAS2002(Tri &T, Tri &Tref, double (*grad)[3])
{
    const double a3 = -2, a4 = 8;
    const double b1 = 0.7, b2 = 0.75;

    // Tangents
    double a0[3], a1[3];
    for (int j = 0; j < 3; j++) {
        a0[j] = T.vert[0]->x[j] - T.vert[2]->x[j];
        a1[j] = T.vert[1]->x[j] - T.vert[2]->x[j];
    }

    // Invariants
    double A = m_ddot(4, &T.a[0][0], &Tref.arcp[0][0]);
    double B = T.detA/Tref.detA;

    double alpha = sqrt(B) - 1;
    double beta = 0.5*A/sqrt(B) - 1;

    if (grad) {
	double dA[2][3] = {0.0}, dB[2][3] = {0.0};
	for (int i = 0; i < 2; i++) {
	    m_daxpy(3, 2*Tref.arcp[i][0], a0, dA[i]);
	    m_daxpy(3, 2*Tref.arcp[i][1], a1, dA[i]);

	    m_daxpy(3, 2*T.detA/Tref.detA*T.arcp[i][0], a0, dB[i]);
	    m_daxpy(3, 2*T.detA/Tref.detA*T.arcp[i][1], a1, dB[i]);
	}

	double dalpha[2][3] = {0.0}, dbeta[2][3] = {0.0};
	for (int i = 0; i < 2; i++) {
	    m_daxpy(3, 0.5/sqrt(B), dB[i], dalpha[i]);

	    m_daxpy(3, 0.5/sqrt(B), dA[i], dbeta[i]);
	    m_daxpy(3, -0.25*A/(B*sqrt(B)), dB[i], dbeta[i]);
	}

	m_dclear(9, *grad);
	for (int i = 0; i < 2; i++) {
	    m_daxpy(3, 
	    	ED*(2*alpha + 3*a3*alpha*alpha + 4*a4*pow(alpha,3))
		+ ES*b1*beta,
		dalpha[i], grad[i]);

	    m_daxpy(3, 
	    	ES*(1 + b1*alpha + 2*b2*beta), 
		dbeta[i], grad[i]);
	}

	for (int j = 0; j < 3; j++) {
	    grad[2][j] = -grad[0][j] - grad[1][j];
	}

	m_dscal(9, Tref.area, *grad);
    }

    return Tref.area*( ED*(alpha*alpha + a3*pow(alpha,3) + a4*pow(alpha,4))
	             + ES*(beta + b1*alpha*beta + b2*beta*beta) );
}


double RedCell::elasticEnergy()
{
    double E = 0.0;

    for (int iface = 0; iface < numFaces(); ++iface) {
        Tri &T = faces[iface];
        Tri &T_ref = cellRef->faces[iface];
        E += tri_elasticity(T, T_ref);
    }

    return E;
}


void RedCell::elasticForce(MArray<double,2> &f)
{
    // Init
    int nvert = numVerts();
    f.resize(nvert,3);
    f = 0.0;

    for (int iface = 0; iface < numFaces(); ++iface) {
	Tri &T = faces[iface];
	Tri &T_ref = cellRef->faces[iface];
	double DW[3][3];

	tri_elasticity(T, T_ref, DW);

	for (int l = 0; l < 3; ++l) {
	    int ivert = T.ivert[l];
	    m_daxpy(3, 1.0/vertArea(ivert), DW[l], &f(ivert,0));
	}
    } // iface
}


double RedCell::bendEnergy()
{
    double E = 0.0;
    double DA = 0.0;	// Area difference between the two lipid leaflets

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

	    E += 2*EB*square(q.H - H0)*q.area;
	    DA += 2*(q.H - H0_GLB)*q.area;
	} // iq
    } // iface

    E += EB_GLB*DA*DA/areaTar;	// add the global bending energy
    return E;
}


/* Bending force 
 * Arguments:
 *  f -- the bending force density
 * Note:
 *  Prerequisites: vertArea, qpoints */
void RedCell::bendForce(MArray<double,2> &f)
{
    int nvert = numVerts();
    int nface = numFaces();

    // lhs = lumped lhs stiffness matrix
    MArray<double,1> lhs(nvert);
    lhs = 0.0;

    double DA = 0.0;	// area difference between the two leaflets

    // gW -- gradient of local bending energy
    // gDA -- gradient of DA
    MArray<double,2> gW(nvert,3), gDA(nvert,3);
    gW = 0.0;
    gDA = 0.0;

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

	MArray<double,2> gW_loc(N,3), gDA_loc(N,3);
	gW_loc = 0.0;
	gDA_loc = 0.0;

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

		lhs_loc(alph) += f(alph)*q.area;

		m_daxpy(3, 2*2*(q.H - H0)*q.area, d_H, &gW_loc(alph,0));
		m_daxpy(3, 2*square(q.H - H0)*q.area, d_J, &gW_loc(alph,0));

		m_daxpy(3, 2*q.area, d_H, &gDA_loc(alph,0));
		m_daxpy(3, 2*(q.H-H0_GLB)*q.area, d_J, &gDA_loc(alph,0));
	    } // alph

	    DA += 2*(q.H - H0_GLB)*q.area;
	} // iq

	// Add to global lhs and rhs
	for (int alph = 0; alph < N; alph++) {
	    int i = ctrlPts[iface][alph];

	    lhs(i) += lhs_loc(alph);

	    m_dadd(3, &gW_loc(alph,0), &gW(i,0));
	    m_dadd(3, &gDA_loc(alph,0), &gDA(i,0));
	}
    } // iface

    // Convert gradient of surface energy to force density
    f.resize(nvert,3);
    f = 0.0;

    for (int ivert = 0; ivert < nvert; ivert++) {
        m_daxpy(3, EB, &gW(ivert,0), &f(ivert,0));
	m_daxpy(3, EB_GLB*2*DA/areaTar, &gDA(ivert,0), &f(ivert,0));

	m_dscal(3, 1.0/vertArea(ivert), &f(ivert,0));
    }
} 


/*
void RedCell::tensionForce(MArray<double,2> &f)
{
    int nvert = numVerts();
    int nface = numFaces();

    f.resize(nvert,3);
    f = 0.0;

    // Tension coefficient
    double sigma = ED_GLB/areaTar*(area/areaTar - 1);

    for (int iface = 0; iface < nface; iface++) {
	Tri &face = faces[iface];
	double dA[3][3];
	face.gradArea(dA);

	for (int l = 0; l < 3; l++) {
	   int ivert = face.ivert[l];
	   FOR_J3 f(ivert,j) += sigma*dA[l][j];
	}
    } // iface

    for (int ivert = 0; ivert < nvert; ivert++) {
        m_dscal(3, 1.0/vertArea(ivert), &f(ivert,0));
    }
} */


/* Surface tension force
 * Arguments:
 *   sigma -- the surface tension
 *   f -- the surface tension force */
void RedCell::tensionForce(double sigma, MArray<double,2> &f)
{
    int nvert = numVerts();
    int nface = numFaces();

    f.resize(nvert,3);
    f = 0.0;

    for (int iface = 0; iface < nface; iface++) {
	Tri &face = faces[iface];
	double dA[3][3];
	face.gradArea(dA);

	for (int l = 0; l < 3; l++) {
	   int ivert = face.ivert[l];
	   FOR_J3 f(ivert,j) += sigma*dA[l][j];
	}
    } // iface

    for (int ivert = 0; ivert < nvert; ivert++) {
        m_dscal(3, 1.0/vertArea(ivert), &f(ivert,0));
    }
} 


/* Membrane viscous force 
 * Arguments:
 *   v -- the surface velocity 
 *   mu -- the membrane viscosity
 *   f -- the density of the viscous force */
/*
void RedCell::viscousForce(MArray<double,2> &v, double mu, MArray<double,2> &f)
{
    int nvert = numVerts();
    f.resize(nvert,3);
    f = 0.0;

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &face = faces[iface];

        double vtmp[3][3];
	for (int l = 0; l < 3; l++) {
	    int ivert = face.ivert[l];
	    m_dcopy(3, &v(ivert,0), vtmp[l]);

	    double vn = m_ddot(3, vtmp[l], face.normal);
	    m_daxpy(3, -vn, face.normal, vtmp[l]);
        }

	double gradv[2][3];
	FOR_D3 {
	    gradv[0][d] = vtmp[0][d] - vtmp[2][d];
	    gradv[1][d] = vtmp[1][d] - vtmp[2][d];
	}

	double a0[3], a1[3], a0_rcp[3], a1_rcp[3];
	FOR_D3 {
	    a0[d] = face.vert[0]->x[d] - face.vert[2]->x[d];
	    a1[d] = face.vert[1]->x[d] - face.vert[2]->x[d];
        }

	FOR_D3 {
	    a0_rcp[d] = face.arcp[0][0]*a0[d] + face.arcp[0][1]*a1[d];
	    a1_rcp[d] = face.arcp[1][0]*a0[d] + face.arcp[1][1]*a1[d];
	}

	// div = in-plane divergence
	double div = m_ddot(3, gradv[0], a0_rcp)
	           + m_ddot(3, gradv[1], a1_rcp);

	// tau = in-plane viscous stress
	double tau[3][3] = {0.0};
	for (int i = 0; i < 3; i++) 
	for (int j = 0; j < 3; j++) {
	    tau[i][j] += gradv[0][i]*a0_rcp[j] + gradv[1][i]*a1_rcp[j];
	    tau[i][j] += gradv[0][j]*a0_rcp[i] + gradv[1][j]*a1_rcp[i];
	    if (i == j) tau[i][j] -= div;
	}

	for (int l = 0; l < 3; l++) {
            int iv0 = face.ivert[l];
            int iv1 = face.ivert[(l+1)%3];
            int iv2 = face.ivert[(l+2)%3];

            double t[3];
            FOR_D3 t[d] = 0.5*(verts[iv2].x[d] - verts[iv1].x[d]);

            double txn[3];
            cross_product(t, face.normal, txn);

            FOR_I3 FOR_J3 f(iv0,i) += tau[i][j]*txn[j];
        }
    }

    for (int ivert = 0; ivert < nvert; ivert++) {
        m_dscal(3, 1.0/vertArea(ivert), &f(ivert,0));
    }

    f *= mu;
} */


/* Calculate the membrane viscous force
 * Arguments:
 *   v -- the surface velocity
 *   f -- the redisual viscous stress
 * On return:
 *   The dissipation rate */
double RedCell::viscousForce(MArray<double,2> &v, MArray<double,2> &f)
{
    // Init
    int nvert = numVerts();
    f.resize(nvert,3);
    f = 0.0;

    double dsp = 0.0;

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &face = faces[iface];

	// Define the two tangential directions
	double a0[3], a1[3];
	FOR_D3 a0[d] = face.vert[1]->x[d] - face.vert[0]->x[d];
	normalizeVec3D(a0);
	cross_product(face.normal, a0, a1);

	// Map things to the tangent plane
        double xtmp[3][2], vtmp[3][2];
	for (int l = 0; l < 3; l++) {
	    int ivert = face.ivert[l];

	    xtmp[l][0] = m_ddot(3, a0, face.vert[l]->x);
	    xtmp[l][1] = m_ddot(3, a1, face.vert[l]->x);

	    vtmp[l][0] = m_ddot(3, a0, &v(ivert,0));
	    vtmp[l][1] = m_ddot(3, a1, &v(ivert,0));
        }

	double gradv[2][2], lhs[2][2], rhs[2][2];

	for (int i = 0; i < 2; i++) {
	    lhs[i][0] = xtmp[0][i] - xtmp[2][i];
	    lhs[i][1] = xtmp[1][i] - xtmp[2][i];

	    rhs[i][0] = vtmp[0][i] - vtmp[2][i];
	    rhs[i][1] = vtmp[1][i] - vtmp[2][i];
	}

	double ilhs[2][2] = {0.0};
	invMat2(lhs, ilhs);

	for (int i = 0; i < 2; i++)
	for (int j = 0; j < 2; j++) {
	    gradv[i][j] = rhs[i][0]*ilhs[0][j] + rhs[i][1]*ilhs[1][j];
	}

	double div = gradv[0][0] + gradv[1][1];

	double tauTmp[2][2];
	for (int i = 0; i < 2; i++)
	for (int j = 0; j < 2; j++) {
	    tauTmp[i][j] = gradv[i][j] + gradv[j][i];
	    if (i == j) tauTmp[i][j] -= div;
	}

	// Add dissipation
	dsp += face.area*m_ddot(4, *tauTmp, *gradv);

	// Add the contribution to the viscous residual force
	// at each vertex
	for (int l = 0; l < 3; l++) {
            double t[2];
	    for (int i = 0; i < 2; i++) {
                t[i] = 0.5*(xtmp[(l+2)%3][i] - xtmp[(l+1)%3][i]);
            }

            double txn[2] = { t[1], -t[0] };

	    double ftmp[2];
	    for (int i = 0; i < 2; i++) {
	        ftmp[i] = tauTmp[i][0]*txn[0] + tauTmp[i][1]*txn[1];
	    }

	    // Map back to 3D
	    int ivert = face.ivert[l];
            FOR_I3 f(ivert,i) += ftmp[0]*a0[i] + ftmp[1]*a1[i];
        }
    }

    for (int ivert = 0; ivert < nvert; ivert++) {
        m_dscal(3, 1.0/vertArea(ivert), &f(ivert,0));
    }

    // What we are calculating is the reaction force from the fluid to
    // the membrane
    f *= -muM;
    dsp *= muM;

    return dsp;
}


/* Calculate the norm of the strain rate
 * Arguments:
 *  v -- the velocity
 *  sl2 -- L2 norm of the strain rate */
void RedCell::strainRate(MArray<double,2> &v, MArray<double,1> &s2)
{
    int nvert = numVerts();
    s2.resize(nvert);
    s2 = 0.0;

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &face = faces[iface];

	// Define the two tangential directions
	double a0[3], a1[3];
	FOR_D3 a0[d] = face.vert[1]->x[d] - face.vert[0]->x[d];
	normalizeVec3D(a0);
	cross_product(face.normal, a0, a1);

	// Map things to the tangent plane
        double xtmp[3][2], vtmp[3][2];
	for (int l = 0; l < 3; l++) {
	    int ivert = face.ivert[l];

	    xtmp[l][0] = m_ddot(3, a0, face.vert[l]->x);
	    xtmp[l][1] = m_ddot(3, a1, face.vert[l]->x);

	    vtmp[l][0] = m_ddot(3, a0, &v(ivert,0));
	    vtmp[l][1] = m_ddot(3, a1, &v(ivert,0));
        }

	double gradv[2][2], lhs[2][2], rhs[2][2];

	for (int i = 0; i < 2; i++) {
	    lhs[i][0] = xtmp[0][i] - xtmp[2][i];
	    lhs[i][1] = xtmp[1][i] - xtmp[2][i];

	    rhs[i][0] = vtmp[0][i] - vtmp[2][i];
	    rhs[i][1] = vtmp[1][i] - vtmp[2][i];
	}

	double ilhs[2][2] = {0.0};
	invMat2(lhs, ilhs);

	for (int i = 0; i < 2; i++)
	for (int j = 0; j < 2; j++) {
	    gradv[i][j] = rhs[i][0]*ilhs[0][j] + rhs[i][1]*ilhs[1][j];
	}

	gradv[0][1] = gradv[1][0] = 0.5*(gradv[0][1] + gradv[1][0]);
	double div = gradv[0][0] + gradv[1][1];
	for (int i = 0; i < 2; i++) gradv[i][i] -= 0.5*div;

	double s2tmp = face.area*m_ddot(4, *gradv, *gradv);

	for (int l = 0; l < 3; l++) {
	    int ivert = face.ivert[l];
	    s2(ivert) += THRD*s2tmp;
        }
    }

    for (int ivert = 0; ivert < nvert; ivert++) {
        s2(ivert) /= vertArea(ivert);
    }
}


/* Tweezer force 
 * Argument:
 *  Ftot -- total force
 *  f -- tweezer force density */
void RedCell::tweezerForce(double Ftot, MArray<double,2> &f)
{
    // Contact area has a diameter of 2 um (Mills and Suresh, MCB, 2004)
    // Length scale is 2.82 um
    // Therefore, 10.7% among all the mesh points are in the contact region    
    int nvert = numVerts();
    int N = (int)round(0.5*0.107*nvert);

    // Apply the tweezer force on N mesh points on each end
    f.resize(nvert,3);
    f = 0.0;

    vector<double> x(nvert);
    for (int ivert = 0; ivert < nvert; ivert++) {
        x[ivert] = verts[ivert].x[0];
    }
    sort( x.begin(), x.end() );

    double x0 = x[N];
    double x1 = x[nvert-N];
    if (x0 > x1) swap(x0, x1);

    // contact = -1 for x < x0
    //         = 0  if no contact
    //         = 1  for x > x1
    MArray<int,1> contact(nvert);
    contact = 0;

    double A0 = 0.0;
    double A1 = 0.0;
    for (int ivert = 0; ivert < nvert; ivert++) {
	if (verts[ivert].x[0] < x0) {
	    contact(ivert) = -1;
	    A0 += vertArea(ivert);
	} else if (verts[ivert].x[0] > x1) {
	    contact(ivert) = 1;
	    A1 += vertArea(ivert);
	}
    }

    // We are computing the hydrodynamic force on the surface,
    // so it is in the opposite direction of the external
    // tweezing forces
    for (int ivert = 0; ivert < nvert; ivert++) {
        if (contact(ivert) == -1) {
	    f(ivert,0) = Ftot/A0;
	} else if (contact(ivert) == 1) {
	    f(ivert,0) = -Ftot/A1;
	}
    }
}


/* Differentiate the elastic and bend force
 * Arguments:
 *   dx -- the displacement
 *   df -- (df/dx)*dx, the change in elastic and bend force */
void RedCell::diffElasBendForce(MArray<double,2> &dx, MArray<double,2> &df)
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

    // Unperturbed force
    MArray<double,2> felas_old(nvert,3), fbend_old(nvert,3), ftens_old;
    elasticForce(felas_old);
    bendForce(fbend_old);
//    tensionForce(ftens_old);

    // Perturbed force
    MArray<double,2> felas_new(nvert,3), fbend_new(nvert,3), ftens_new;

    for (int ivert = 0; ivert < nvert; ivert++) {
        Point &vert = verts[ivert];
	m_daxpy(3, scal, &dx(ivert,0), vert.x);
    }
    updateGeometry();

    elasticForce(felas_new);
    bendForce(fbend_new);
//    tensionForce(ftens_new);

    // Restore shape
    setCoords(xsave);
    updateGeometry();

    // Use finite-difference
    df.resize(nvert,3);

    for (int ivert = 0; ivert < nvert; ivert++) {
        FOR_J3 df(ivert,j) = 
	    (felas_new(ivert,j) - felas_old(ivert,j))
	  + (fbend_new(ivert,j) - fbend_old(ivert,j));
//	  + (ftens_new(ivert,j) - ftens_old(ivert,j));
    }
    df *= 1.0/scal;
}


/* Rebuild vertex valences and build control point list
 * Note:
 *  -- This needs to be done once and only once when the mesh
 *     topology is changed */
void RedCell::buildControlPoints()
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


void RedCell::calcQuadPoint(
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
