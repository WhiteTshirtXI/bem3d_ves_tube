#include "cxxheaders.h"
#include "greenfunc.h"
#include "tri.h"
#include "mathfunc.h"

/* void greenf_unbound 
 * void greenf_plate 
 * void greenf_int_on_tri_regular 
 * void greenf_int_on_tri_singular */


/* Calculate the velocity induced by a stokeslet in an infinite space
 * Arguments:
 *   x0 -- 
 *   x -- 
 *   n -- normal direction at x
 *   G -- Stokeslet
 *   f -- 
 *   fG -- fG_j = f_i(x) * G_ij(x,x0)
 *   T -- T(x,x0)_ij = T(x,x0)_ijk n_k(x), Stresslet*normal
 *   g --
 *   gT -- gT_j = g_i(x) T_ijk(x,x0) n_k(x) */
void greenf_unbound(const double *x0, const double *x, const double *n,
		double (*G)[3], const double *f, double *fG, 
		double (*T)[3], const double *g, double *gT)
{
    // parameter check and initialize
    assert(x0);
    assert(x);

    if (G) m_dclear(9, *G);
    if (fG) {
        assert(f);
        m_dclear(3, fG);
    }

    if (T) {
        assert(n);
	m_dclear(9, *T);
    }
    if (gT) {
        assert(g);
	assert(n);
	m_dclear(3, gT);
    }

    double xx[3], r, r2, ir, ir2, ir3, ir5;

    FOR_I3 xx[i] = x[i] - x0[i];
    r2 = m_ddot(3, xx, xx);
    r = sqrt(r2);
    ir = 1.0/r;
    ir2 = 1.0/r2;
    ir3 = ir*ir2;
    ir5 = ir2*ir3;

    if (G) {
        FOR_I3 {
	    G[i][i] += ir;
	    FOR_J3 G[i][j] += ir3*xx[i]*xx[j];
	}
    }

    if (fG) {
        double xxf = m_ddot(3, xx, f);
        FOR_J3 fG[j] += ir*f[j] + ir3*xxf*xx[j];
    }

    if (T) {
        double xxn = m_ddot(3,xx,n);
	double c = -6*ir5*xxn;
	FOR_I3 {
	    FOR_J3 T[i][j] += c*xx[i]*xx[j];
	}
    }

    if (gT) {
        double gxx = m_ddot(3, g, xx);
        double xxn = m_ddot(3,xx,n);
	double c = -6*ir5*gxx*xxn;
	FOR_J3 gT[j] += c*xx[j];
    }
}


/* Green's function for wall bounded flow
 * Note:
 *   1. The wall is the z=0 plane
 * Reference:
 *   1. C. Pozrikids, "Boundary integral and isngularity methods for 
 *   linearized viscous flow", Cambridge, 1992 */
void greenf_plate(const double *x0, const double *x, const double *n,
		double (*G)[3], const double *f, double *fG, 
		double (*T)[3], const double *g, double *gT)
{
    // check parameters and initialize
    assert(x0);
    assert(x);

    if (G) m_dclear(9, *G);
    if (fG) {
        assert(f);
        m_dclear(3, fG);
    }

    if (T) {
        assert(n);
	m_dclear(9, *T);
    }
    if (gT) {
        assert(g); 
	assert(n);
	m_dclear(3, gT);
    }

    double xx[3], r, r2, ir, ir2, ir3, ir5;

    // This part is the same as the unbounded space Green's function
    FOR_I3 xx[i] = x[i] - x0[i];
    r2 = m_ddot(3, xx, xx);
    r = sqrt(r2);
    ir = 1.0/r;
    ir2 = 1.0/r2;
    ir3 = ir*ir2;
    ir5 = ir2*ir3;

    if (G) {
        FOR_I3 {
	    G[i][i] += ir;
	    FOR_J3 G[i][j] += ir3*xx[i]*xx[j];
        }
    }

    if (fG) {
        double xxf = m_ddot(3, xx, f);
        FOR_I3 fG[i] += ir*f[i] + ir3*xx[i]*xxf;
    }

    if (T) {
        double xxn = m_ddot(3,xx,n);
	double c = -6*ir5*xxn;
	FOR_I3 {
	    FOR_J3 T[i][j] += c*xx[i]*xx[j];
	}
    }

    if (gT) {
        double gxx = m_ddot(3, g, xx);
        double xxn = m_ddot(3,xx,n);
	double c = -6*ir5*gxx*xxn;
	FOR_J3 gT[j] += c*xx[j];
    }


    // Contribution from image point
    const double xim[3] = {x0[0], x0[1], -x0[2]};

    FOR_I3 xx[i] = x[i] - xim[i];
    r2 = m_ddot(3, xx, xx);
    r = sqrt(r2);
    ir = 1.0/r;
    ir2 = 1.0/r2;
    ir3 = ir*ir2;
    ir5 = ir2*ir3;

    // -f at image point
    if (G) {
        FOR_I3 {
	    G[i][i] -= ir;
	    FOR_J3 G[i][j] -= ir3*xx[i]*xx[j];
        }
    }

    if (fG) {
        double xxf = m_ddot(3, xx, f);
        FOR_I3 fG[i] -= (ir*f[i] + ir3*xx[i]*xxf);
    }

    if (T) {
        double xxn = m_ddot(3,xx,n);
	double c = -6*ir5*xxn;
	FOR_I3 {
	    FOR_J3 T[i][j] -= c*xx[i]*xx[j];
	}
    }


    if (gT) {
        double gxx = m_ddot(3, g, xx);
        double xxn = m_ddot(3,xx,n);
	double c = -6*ir5*gxx*xxn;
	FOR_J3 gT[j] -= c*xx[j];
    }

    // Now need to flip the normal force component at the image point
    const double h = x[2];
    const double h0 = x0[2];
    const double hh0_3 = h*h0*ir3;	// h*h0/r^3
    const double hh0_5 = hh0_3*ir2;	// h*h0/r^5
    const double hh0_7 = hh0_5*ir2;	// h*h0/r^7
    const double h0_3 = h0*ir3;		// h0/r^3
    const double h0_5 = h0_3*ir2;	// h0/r^5

    if (G) {
        double Gim[3][3];
	m_dclear(9, *Gim);

	FOR_I3 
	FOR_J3 {
	    if (i==j) Gim[i][j] += -2*hh0_3;
	    Gim[i][j] += 6*hh0_5*xx[i]*xx[j];
	    if (j==2) Gim[i][j] += -2*h0_3*xx[i];
	    if (i==2) Gim[i][j] += 2*h0_3*xx[j];
	}

	FOR_I3 {
	    G[i][0] += Gim[i][0];
	    G[i][1] += Gim[i][1];
	    G[i][2] -= Gim[i][2];	// effect of flip
	}
    }

    if (fG) {
	double fGim[3] = {0.0, 0.0, 0.0};
        double xxf = m_ddot(3, xx, f);

        FOR_J3 {
	    fGim[j] += -2*hh0_3*f[j] + 6*hh0_5*xxf*xx[j];
	    if (j==2) fGim[j] += -2*h0_3*xxf;
	    fGim[j] += 2*h0_3*f[2]*xx[j];
	}

	fG[0] += fGim[0];
	fG[1] += fGim[1];
	fG[2] -= fGim[2];
    }

    if (T) {
        double Tim[3][3];
	m_dclear(9, *Tim);
	double xxn = m_ddot(3, xx, n);

        FOR_I3
	FOR_J3 {
	    Tim[i][j] += 12*hh0_5*(n[i]*xx[j] + xx[i]*n[j]);
	    if (i==j) Tim[i][j] += 12*hh0_5*xxn;
	    Tim[i][j] -= 60*hh0_7*xx[i]*xx[j]*xxn;

	    Tim[i][j] -= 12*h0_5*n[i]*xx[j]*xx[2];
	    if (j==2) Tim[i][j] += 12*h0_5*xx[i]*xxn;
	}

	FOR_I3 {
	    T[i][0] += Tim[i][0];
	    T[i][1] += Tim[i][1];
	    T[i][2] -= Tim[i][2];	// flip
	}
    }


    if (gT) {
        double gTim[3] = {0.0, 0.0, 0.0};
	double xxg = m_ddot(3, xx, g);
	double xxn = m_ddot(3, xx, n);
	double gn = m_ddot(3, g, n);

	FOR_J3 {
	    gTim[j] += 12*hh0_5*(gn*xx[j] + g[j]*xxn + xxg*n[j]);
	    gTim[j] -= 60*hh0_7*(xxg*xx[j]*xxn);
	    gTim[j] -= 12*h0_5*gn*xx[j]*xx[2];
	    if (j==2) gTim[j] += 12*h0_5*xxg*xxn;
	}

	gT[0] += gTim[0];
	gT[1] += gTim[1];
	gT[2] -= gTim[2];
    }
}


/* Regular surface integral of a Green's function kernel over a triangle
 * Arguments:
 *   x0 -- the target point
 *   tri -- the triangle
 *   greenf --
 *   c1 -- single layer coefficient
 *   rhs1 -- 
 *   lhs1 --
 *   c2 -- double layer coefficient
 *   rhs2 --
 *   lsh2 -- 
 *   qrule -- qudrature rule
 *   Note: 
 *    1. Assume the single-layer density is tri->vert[l]->f, l=0,1,2
 *    2. Assume the double-layer density is tri->vert[l]->g, l=0,1,2
 *    3. For single layer
 *        int f_i(x) G_ij(x,x0) dS(x) 
 *          = sum_l lhs[j][l][i] * f[l][i]
 *          = rhs[j]
 *       where l=0,1,2 is the vertex index
 *    4. For double layer
 *        int g_i(x) T_ijk(x,x0) n_k(x) dS(x)
 *          = sum_l lhs[j][l][i] * g[l][i]
 *          = rhs[j] */
void greenf_int_on_tri_regular(pGreenFunc greenf, const double *x0, const Tri &tri,
		double c1, double *rhs1, double (*lhs1)[3][3],
		double c2, double *rhs2, double (*lhs2)[3][3],
		const char *qrule)
{
    // initialize
    if (rhs1) m_dclear(3, rhs1);
    if (lhs1) m_dclear(27, **lhs1);

    if (rhs2) m_dclear(3, rhs2);
    if (lhs2) m_dclear(27, **lhs2);

    Quad2D &q = quadrature::select_rule_2d(qrule);

    // Copy variable
    double xtri[3][3], ftri[3][3], gtri[3][3];

    FOR_I3 m_dcopy(3, tri.vert[i]->x, xtri[i]);

    if (rhs1) 
        FOR_I3 m_dcopy(3, tri.vert[i]->f, ftri[i]);
    else
        m_dclear(9, *ftri);

    if (rhs2)
        FOR_I3 m_dcopy(3, tri.vert[i]->g, gtri[i]);
    else
        m_dclear(9, *gtri);

    for (int iq = 0; iq < q.n(); ++iq) {
        double s, t, w0, w1, w2;
        double xq[3], fq[3], gq[3], dA;

        s = q.x(iq);
        t = q.y(iq);

        w0 = 1 - s - t;
        w1 = s;
        w2 = t;

        dA = tri.detJ*q.w(iq);

	FOR_I3 {
	    xq[i] = w0*xtri[0][i] + w1*xtri[1][i] + w2*xtri[2][i];

	    fq[i] = w0*ftri[0][i] + w1*ftri[1][i] + w2*ftri[2][i];
	    fq[i] *= dA;

	    gq[i] = w0*gtri[0][i] + w1*gtri[1][i] + w2*gtri[2][i];
	    gq[i] *= dA;
        }

	double fG[3], G[3][3], *pfG, (*pG)[3];
	double gT[3], T[3][3], *pgT, (*pT)[3];

	pfG = rhs1?  fG : NULL;
	pG = lhs1? G : NULL;

	pgT = rhs2?  gT : NULL;
	pT = lhs2? T : NULL;

	(*greenf)(x0, xq, tri.normal, pG, fq, pfG, pT, gq, pgT);

	if (rhs1) m_daxpy(3, c1, fG, rhs1);
	if (rhs2) m_daxpy(3, c2, gT, rhs2);

	if (lhs1) {
	    FOR_I3 {
	        FOR_J3 {
	            lhs1[i][0][j] += c1*w0*G[j][i]*dA;
	            lhs1[i][1][j] += c1*w1*G[j][i]*dA;
	            lhs1[i][2][j] += c1*w2*G[j][i]*dA;
	        } 
	    } 
	}

	if (lhs2) {
	    FOR_I3 {
	        FOR_J3 {
	            lhs2[i][0][j] += c2*w0*T[j][i]*dA;
	            lhs2[i][1][j] += c2*w1*T[j][i]*dA;
	            lhs2[i][2][j] += c2*w2*T[j][i]*dA;
	        } 
	    } 
	}
    } //iq
}


/* Calculate the singular surface integral
 * Arguments:
 *   All but nqp -- Same arguments as those in greenf_int_on_tri_regular
 *   nq -- nubmer of qudrature points for Duffy integration in each direction */
void greenf_int_on_tri_singular(pGreenFunc greenf, const double *x0, const Tri &tri, 
		double c1, double *rhs1, double (*lhs1)[3][3],
		double c2, double *rhs2, double (*lhs2)[3][3],
		int nq)
{
    // initialize
    if (rhs1) m_dclear(3, rhs1);	
    if (lhs1) m_dclear(27, **lhs1);
    if (rhs2) m_dclear(3, rhs2);	
    if (lhs2) m_dclear(27, **lhs2);

    double xtri[3][3], ftri[3][3], gtri[3][3];

    FOR_I3 m_dcopy(3, tri.vert[i]->x, xtri[i]);

    if (rhs1) 
        FOR_I3 m_dcopy(3, tri.vert[i]->f, ftri[i]);
    else
        m_dclear(9, *ftri);

    if (rhs2)
        FOR_I3 m_dcopy(3, tri.vert[i]->g, gtri[i]);
    else
        m_dclear(9, *gtri);

    // Find projection
    double s0, t0;
    double dist0 = minDistToTri(x0, xtri, s0, t0);

    // Set up Duffy quadrature rule
    double *sq = new double[nq]; 
    double *wq = new double[nq];
    gauleg(0.0, 1.0, nq, sq, wq);

    // Divide the triangle into 3 sub-triangles
    for (int isub = 0; isub < 3; ++isub) {
        double detjSub, s1, s2, t1, t2;

        // (s0,t0), (s1,t1), (s2,t2) are the reference coordinates of
        //  the three corners of the sub-triangle in the original triangle
        switch (isub) {
            case 0:
                s1 = 1.0;    t1 = 0.0;
                s2 = 0.0;    t2 = 1.0;
                detjSub = (1 - s0 - t0)*tri.detJ;
                break;
            case 1:
                s1 = 0.0;    t1 = 1.0;
                s2 = 0.0;    t2 = 0.0;
                detjSub = s0*tri.detJ;
                break;
            case 2:
                s1 = 0.0;    t1 = 0.0;
                s2 = 1.0;    t2 = 0.0;
                detjSub = t0*tri.detJ;
                break;
        }

        if (fabs(detjSub) < 1.E-10*tri.detJ) continue;

        for (int iq = 0; iq < nq; ++iq)
        for (int jq = 0; jq < nq; ++jq) {
            double sloc, tloc, s, t, xq[3], fq[3], gq[3], dA;
            double w0, w1, w2;

            sloc = sq[iq];
            tloc = sloc*sq[jq];
            dA = detjSub*sloc*wq[iq]*wq[jq];

            s = (1 - sloc)*s0 + (sloc - tloc)*s1 + tloc*s2;
            t = (1 - sloc)*t0 + (sloc - tloc)*t1 + tloc*t2;

            w0 = 1 - s - t;
            w1 = s;
            w2 = t;

	    FOR_I3 {
	        xq[i] = w0*xtri[0][i] + w1*xtri[1][i] + w2*xtri[2][i];

		fq[i] = w0*ftri[0][i] + w1*ftri[1][i] + w2*ftri[2][i];
		fq[i] *= dA;

		gq[i] = w0*gtri[0][i] + w1*gtri[1][i] + w2*gtri[2][i];
		gq[i] *= dA;
	    }

	    double fG[3], G[3][3], *pfG, (*pG)[3];
	    double gT[3], T[3][3], *pgT, (*pT)[3];

	    pfG = rhs1?  fG : NULL;
	    pG = lhs1? G : NULL;

	    pgT = rhs2?  gT : NULL;
	    pT = lhs2? T : NULL;

	    (*greenf)(x0, xq, tri.normal, pG, fq, pfG, pT, gq, pgT);

	    if (rhs1) m_daxpy(3, c1, pfG, rhs1);
	    if (rhs2) m_daxpy(3, c2, pgT, rhs2);

	    if (lhs1) {
		FOR_I3 {
		    FOR_J3 {
			lhs1[i][0][j] += c1*w0*G[j][i]*dA;
			lhs1[i][1][j] += c1*w1*G[j][i]*dA;
			lhs1[i][2][j] += c1*w2*G[j][i]*dA;
		    } 
		} 
	    }

	    if (lhs2) {
		FOR_I3 {
		    FOR_J3 {
			lhs2[i][0][j] += c2*w0*T[j][i]*dA;
			lhs2[i][1][j] += c2*w1*T[j][i]*dA;
			lhs2[i][2][j] += c2*w2*T[j][i]*dA;
		    } 
		} 
	    }

        } //iq,jq
    } //isub

    delete [] sq;
    delete [] wq;
} 
