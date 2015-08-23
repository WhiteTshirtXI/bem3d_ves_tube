#include "cxxheaders.h"
#include "ewald.h"
#include "mathfunc.h"
#include "mesh.h"
#include "debugfunc.h"
#include "quadrature.h"

/* ewald::phys::coeff_exact
 * ewald::phys::coeff
 * ewald::phys::int_on_tri_regular
 * ewald::phys::int_on_tri_singular
 * ewald::phys::int_on_tri_nearsingular
 * ewald::phys::addSurfInt
 * ewald::phys::calcSurfIntMat */


/* Compute Ewald physical sum coefficients exactly
 * Arguments:
 *  r -- separation distance
 *  S1,S2 -- single-layer coefficient
 *  D1 -- double-layer coefficient
 * Note:
 *   -- For single layer, the velocity induced is
 *         S1*x*(x*f) + S2*f, where x is the UNIT vector from source to target
 *   -- For double layer, the velocity induced is
 *         D*x*x*(x*g) */
void ewald::phys::coeff_exact(double r, double *S1, double *S2, double *D1)
{
    // Free space Green's function
    if (alpha <= 0) {
        double ir = 1.0/r;
        if (S1) *S1 = ir;
        if (S2) *S2 = ir;
        if (D1) *D1 = -6.0*ir*ir;
        return;
    }

    double r2, ir, ir2, r2_t, r_t;
    r2 = r*r;
    ir = 1.0/r;
    ir2 = ir*ir;
    r2_t = (M_PI/alpha)*r2;
    r_t = sqrt(r2_t);

    if (S1 || S2) {
        double c1 = erfc(r_t);
        double c2 = 2.0/sqrt(alpha)*exp(-r2_t)*r;

        if (S1) *S1 = (c1 + c2)*ir;
        if (S2) *S2 = (c1 - c2)*ir;
    }

    if (D1) {
        *D1 = -exp(-r2_t)*(12*r_t + 8*r_t*r2_t)/sqrt(M_PI) - 6*erfc(r_t);
        *D1 *= ir2;
    }
}


/* Compute Ewald physical sum coefficients using lookup tables with
 * linear interpolation
 * Arguments:
 *   The same as coeff_exact */
void ewald::phys::coeff(double r, double *S1, double *S2, double *D1)
{
    const int N_tab = 1024;
    static double coeff_tab[N_tab+1][3];
    static double rmax, dr, idr;
    static bool inited = false;

    if (alpha > 1.E10 || alpha <= 0) {
        coeff_exact(r, S1, S2, D1);
        return;
    }

    // Build table
    if (! inited ) {
        rmax = 2.5*sqrt(alpha);
        dr = rmax/N_tab;
        idr = 1.0/dr;

        // Use limiting values for r = 0
        coeff_tab[0][0] = 1.0;
        coeff_tab[0][1] = 1.0;
        coeff_tab[0][2] = -6.0;

        for (int i = 1; i <= N_tab; i++) {
            double r = i*dr;
            coeff_exact(r, &coeff_tab[i][0], &coeff_tab[i][1], &coeff_tab[i][2]);

            coeff_tab[i][0] *= r;
            coeff_tab[i][1] *= r;
            coeff_tab[i][2] *= (r*r);
        }

        inited = true;
    }

    if (r >= rmax) {
        if (S1) *S1 = 0.0;
        if (S2) *S2 = 0.0;
        if (D1) *D1 = 0.0;
        return;
    } else {
        double s = r*idr;
        int i = min((int)floor(s), N_tab-1);
        double w0 = i + 1 - s;
        double w1 = 1 - w0;

        double ir = 1.0/r;
        if (S1) {
            *S1 = w0*coeff_tab[i][0] + w1*coeff_tab[i+1][0];
            *S1 *= ir;
        }
        if (S2) {
            *S2 = w0*coeff_tab[i][1] + w1*coeff_tab[i+1][1];
            *S2 *= ir;
        }
        if (D1) {
            *D1 = w0*coeff_tab[i][2] + w1*coeff_tab[i+1][2];
            *D1 *= ir*ir;
        }
    }
}


/* Regular integration on a triangle
 * Arguments:
 *   xtar -- the target point
 *   tri -- the triangular surface
 *   lhs1 -- single-layer lhs
 *   lhs2 -- double-layer lhs
 *   c1 -- single-layer coeff
 *   c2 -- double-layer coeff
 *   rhs -- the rhs
 *   qrule -- the quadrature rule
 * Note:
 *   1. rhs[d] = c1*lhs1(d,j,k)*f(j,k) + c2*lhs2(d,j,k)*g(j,k)
 *   where f(j,k) is the j-th vertex's k-th component of force density.
 *   lhs2 and g are defined similarly.
 *   2. The target point must be moved to close to the triangle */
void ewald::phys::int_on_tri_regular(const double *xtar, const Tri &tri,
        double (*lhs1)[3][3], double (*lhs2)[3][3],
        double c1, double c2, double *rhs,
        const char *qrule)
{
    bool c1_flag = fabs(c1) > 1.E-10;
    bool c2_flag = fabs(c2) > 1.E-10;

    // Copy variables
    double y[3][3], f[3][3] = {0.0}, g[3][3] = {0.0};

    for (int i = 0; i < 3; ++i) {
        m_dcopy(3, tri.vert[i]->x, y[i]);
        if (c1_flag) m_dcopy(3, tri.vert[i]->f, f[i]);
        if (c2_flag) m_dcopy(3, tri.vert[i]->g, g[i]);
    }

    // Choose 3-point quadrature rule
    Quad2D &q = quadrature::select_rule_2d(qrule);

    // Initialize
    if (lhs1) m_dclear(27, &lhs1[0][0][0]);
    if (lhs2) m_dclear(27, &lhs2[0][0][0]);
    if (rhs) m_dclear(3, rhs);

    for (int iq = 0; iq < q.n(); ++iq) {
        double s, t, w0, w1, w2, dA;
        double xx[3], rr;
        double S1, S2, D1;

        s = q.x(iq);
        t = q.y(iq);

        w0 = 1 - s - t;
        w1 = s;
        w2 = t;

        dA = tri.detJ*q.w(iq);

        for (int ii = 0; ii < 3; ++ii)
            xx[ii] = w0*y[0][ii] + w1*y[1][ii] + w2*y[2][ii] - xtar[ii];

        rr = m_dnrm2(3, xx);
        m_dscal(3, 1.0/rr, xx);
        ewald::phys::coeff(rr, &S1, &S2, &D1);

        if (lhs1) {
            double lhstmp[3][3] = {0.0};
            double fac;

            fac = S1*dA;
            for (int ii = 0; ii < 3; ++ii)
                for (int jj = 0; jj < 3; ++jj) {
                    lhstmp[ii][jj] = fac*xx[ii]*xx[jj];
                }

            fac = S2*dA;
            for (int ii = 0; ii < 3; ++ii) lhstmp[ii][ii] += fac;

            for (int d = 0; d < 3; ++d)
                for (int jj = 0; jj < 3; ++jj) {
                    lhs1[d][0][jj] += w0*lhstmp[d][jj];
                    lhs1[d][1][jj] += w1*lhstmp[d][jj];
                    lhs1[d][2][jj] += w2*lhstmp[d][jj];
                }
        }

        if (lhs2) {
            double lhstmp[3][3] = {0.0};
            double fac;

            fac = D1*m_ddot(3, xx, tri.normal)*dA;

            for (int ii = 0; ii < 3; ++ii)
                for (int jj = 0; jj < 3; ++jj) {
                    lhstmp[ii][jj] = fac*xx[ii]*xx[jj];
                }

            for (int d = 0; d < 3; ++d)
                for (int jj = 0; jj < 3; ++jj) {
                    lhs2[d][0][jj] += w0*lhstmp[d][jj];
                    lhs2[d][1][jj] += w1*lhstmp[d][jj];
                    lhs2[d][2][jj] += w2*lhstmp[d][jj];
                }
        }

        if (rhs) {
            if (c1_flag) {
                // Add single layer potential
                double ff[3], xxff, fac;

                for (int ii = 0; ii < 3; ++ii)
                    ff[ii] = w0*f[0][ii] + w1*f[1][ii] + w2*f[2][ii];

                xxff = m_ddot(3, xx, ff);

                fac = c1*S1*xxff*dA;
                m_daxpy(3, fac, xx, rhs);

                fac = c1*S2*dA;
                m_daxpy(3, fac, ff, rhs);
            }


            if (c2_flag) {
                // Add double layer potential
                double gg[3], xxgg, xxnn, fac;

                for (int ii = 0; ii < 3; ++ii)
                    gg[ii] = w0*g[0][ii] + w1*g[1][ii] + w2*g[2][ii];

                xxgg = m_ddot(3, xx, gg);
                xxnn = m_ddot(3, xx, tri.normal);

                fac = c2*D1*xxnn*xxgg*dA;
                m_daxpy(3, fac, xx, rhs);
            }
        }
    } // iq
}


/* Calculate the integral on a triangle using Duffy rule
 * Arguments:
 *   xtar -- target points
 *   tri -- the triangle
 *   lhs1, lhs2 -- single and double layer matrices
 *   c1, c2 -- single and double layer coefficients
 *   rhs -- the integral
 *   nqp -- the number of quadrature point in Duffy rule
 */
void ewald::phys::int_on_tri_singular(const double *xtar, const Tri &tri,
        double (*lhs1)[3][3], double (*lhs2)[3][3],
        double c1, double c2, double *rhs,
        int nqp)
{
    // Set flag
    bool c1_flag = fabs(c1) > 1.E-10;
    bool c2_flag = fabs(c2) > 1.E-10;

    // Copy variables
    double y[3][3], f[3][3]={}, g[3][3]={};

    for (int i = 0; i < 3; ++i) {
        m_dcopy(3, tri.vert[i]->x, y[i]);
        if (c1_flag) m_dcopy(3, tri.vert[i]->f, f[i]);
        if (c2_flag) m_dcopy(3, tri.vert[i]->g, g[i]);
    }

    // Set up Duffy quadrature rule
    double s0, t0;
    minDistToTri(xtar, y, s0, t0);

    vector<double> sq, tq, wq;
    calcDuffyQuadPoint(s0, t0, nqp, sq, tq, wq);

    // Initialize
    if (lhs1) m_dclear(27, &lhs1[0][0][0]);
    if (lhs2) m_dclear(27, &lhs2[0][0][0]);
    if (rhs) m_dclear(3, rhs);

    for (int iq = 0; iq < sq.size(); iq++) {
	double s, t, xx[3], ff[3], gg[3], rr, dA;
	double w0, w1, w2;
	double S1, S2, D1;

	dA = wq[iq]*tri.detJ;

        s = sq[iq];
	t = tq[iq];

        w0 = 1 - s - t;
        w1 = s;
        w2 = t;

        FOR_I3 xx[i] = w0*y[0][i] + w1*y[1][i] + w2*y[2][i] - xtar[i];
	rr = m_dnrm2(3, xx);
	m_dscal(3, 1.0/rr, xx);
	ewald::phys::coeff(rr, &S1, &S2, &D1);

	if (lhs1) {
	    double lhstmp[3][3] = {0.0};
	    double fac;

	    fac = S1*dA;
	    FOR_I3 FOR_J3 lhstmp[i][j] = fac*xx[i]*xx[j];

	    fac = S2*dA;
	    FOR_I3 lhstmp[i][i] += fac;

	    FOR_I3 
		FOR_J3 {
		    lhs1[i][0][j] += w0*lhstmp[i][j];
		    lhs1[i][1][j] += w1*lhstmp[i][j];
		    lhs1[i][2][j] += w2*lhstmp[i][j];
		}
	} // lhs1

	if (lhs2) {
	    double lhstmp[3][3] = {0.0};
	    double fac;

	    fac = D1*m_ddot(3, xx, tri.normal)*dA;

	    FOR_I3 FOR_J3 lhstmp[i][j] = fac*xx[i]*xx[j];

	    FOR_I3
		FOR_J3 {
		    lhs2[i][0][j] += w0*lhstmp[i][j];
		    lhs2[i][1][j] += w1*lhstmp[i][j];
		    lhs2[i][2][j] += w2*lhstmp[i][j];
		}
	} // lhs2

	if (rhs) {
	    if (c1_flag) {
		double ff[3], xxff, fac;

		FOR_I3 ff[i] = w0*f[0][i] + w1*f[1][i] + w2*f[2][i];
		xxff = m_ddot(3, xx, ff);

		fac = c1*S1*xxff*dA;
		m_daxpy(3, fac, xx, rhs);

		fac = c1*S2*dA;
		m_daxpy(3, fac, ff, rhs);
	    }

	    if (c2_flag) {
		double gg[3], xxnn, xxgg, fac;

		FOR_I3 gg[i] = w0*g[0][i] + w1*g[1][i] + w2*g[2][i];
		xxnn = m_ddot(3, xx, tri.normal);
		xxgg = m_ddot(3, xx, gg);

		fac = c2*D1*xxnn*xxgg*dA;
		m_daxpy(3, fac, xx, rhs);
	    }
	} // rhs
    } // iq
}


/* Calculate the near singular integral on a triangle using log transform
 * Arguments:
 *   xtar -- target points
 *   tri -- the triangle
 *   lhs1, lhs2 -- single and double layer matrices
 *   c1, c2 -- single and double layer coefficients
 *   rhs -- the integral
 *   nqp -- the number of quadrature point in Duffy rule
 */
void ewald::phys::int_on_tri_nearsingular(const double *xtar, const Tri &tri,
        double (*lhs1)[3][3], double (*lhs2)[3][3],
        double c1, double c2, double *rhs,
        int nqp)
{
    // Set flag
    bool c1_flag = fabs(c1) > 1.E-10;
    bool c2_flag = fabs(c2) > 1.E-10;

    // Copy variables
    double y[3][3], f[3][3]={}, g[3][3]={};

    for (int i = 0; i < 3; ++i) {
        m_dcopy(3, tri.vert[i]->x, y[i]);
        if (c1_flag) m_dcopy(3, tri.vert[i]->f, f[i]);
        if (c2_flag) m_dcopy(3, tri.vert[i]->g, g[i]);
    }

    vector<double> sq, tq, dAq;
    calcLogQuadPoint(xtar, tri, nqp, sq, tq, dAq);

    // Initialize
    if (lhs1) m_dclear(27, &lhs1[0][0][0]);
    if (lhs2) m_dclear(27, &lhs2[0][0][0]);
    if (rhs) m_dclear(3, rhs);

    for (int iq = 0; iq < sq.size(); iq++) {
	double s, t, xx[3], ff[3], gg[3], rr, dA;
	double w0, w1, w2;
	double S1, S2, D1;

        s = sq[iq];
	t = tq[iq];
	dA = dAq[iq];

        w0 = 1 - s - t;
        w1 = s;
        w2 = t;

        FOR_I3 xx[i] = w0*y[0][i] + w1*y[1][i] + w2*y[2][i] - xtar[i];
	rr = m_dnrm2(3, xx);
	m_dscal(3, 1.0/rr, xx);
	ewald::phys::coeff(rr, &S1, &S2, &D1);

	if (lhs1) {
	    double lhstmp[3][3] = {0.0};
	    double fac;

	    fac = S1*dA;
	    FOR_I3 FOR_J3 lhstmp[i][j] = fac*xx[i]*xx[j];

	    fac = S2*dA;
	    FOR_I3 lhstmp[i][i] += fac;

	    FOR_I3
		FOR_J3 {
		    lhs1[i][0][j] += w0*lhstmp[i][j];
		    lhs1[i][1][j] += w1*lhstmp[i][j];
		    lhs1[i][2][j] += w2*lhstmp[i][j];
		}
	} // lhs1

	if (lhs2) {
	    double lhstmp[3][3] = {0.0};
	    double fac;

	    fac = D1*m_ddot(3, xx, tri.normal)*dA;

	    FOR_I3 FOR_J3 lhstmp[i][j] = fac*xx[i]*xx[j];

	    FOR_I3
		FOR_J3 {
		    lhs2[i][0][j] += w0*lhstmp[i][j];
		    lhs2[i][1][j] += w1*lhstmp[i][j];
		    lhs2[i][2][j] += w2*lhstmp[i][j];
		}
	} // lhs2

	if (rhs) {
	    if (c1_flag) {
		double ff[3], xxff, fac;

		FOR_I3 ff[i] = w0*f[0][i] + w1*f[1][i] + w2*f[2][i];
		xxff = m_ddot(3, xx, ff);

		fac = c1*S1*xxff*dA;
		m_daxpy(3, fac, xx, rhs);

		fac = c1*S2*dA;
		m_daxpy(3, fac, ff, rhs);
	    }

	    if (c2_flag) {
		double gg[3], xxnn, xxgg, fac;

		FOR_I3 gg[i] = w0*g[0][i] + w1*g[1][i] + w2*g[2][i];
		xxnn = m_ddot(3, xx, tri.normal);
		xxgg = m_ddot(3, xx, gg);

		fac = c2*D1*xxnn*xxgg*dA;
		m_daxpy(3, fac, xx, rhs);
	    }
	} // rhs
    } // iq
}



/* Add surface integration
 * Arguments:
 *   nlist -- the nbrlist
 *   c1, c2 -- the single and double layer coefficients
 *   v -- velocities */
void ewald::phys::addSurfInt(NbrList &nlist, double c1, double c2, double *v)
{
    for (int i = 0; i < nlist.verts.size(); ++i) {
        Point *point = nlist.verts[i];

        for (int p = nlist.firstNbr[i]; p < nlist.firstNbr[i+1]; p++) {
            Tri *T = nlist.faces[p];

            double dist = nlist.dists[p];
            if (dist > rc) continue;

	    // Special case for double-layer potential
	    if (fabs(c1) < 1.E-10 && 
	        ( point == T->vert[0] || 
		  point == T->vert[1] || 
		  point == T->vert[2] ) ) continue;

            // Move target point closer to the triangle
            double xtar[3], vtmp[3];
            m_dcopy(3, point->x, xtar);
            to_CloseBy(T->xc, xtar);

            // Choose integator
            if (dist*dist > T->detJ)
                int_on_tri_regular(xtar, *T, NULL, NULL, c1, c2, vtmp);
            else
                int_on_tri_singular(xtar, *T, NULL, NULL, c1, c2, vtmp);

            m_dadd(3, vtmp, v+3*i);
        } // p
    } // i
}


/* Calculate the lhs matrix due to single- and double- layer operator
 * Arguments:
 *   nlist -- neighbor list
 *   lhs1 -- single-layer matrix
 *   lhs2 -- double-layer matrix */
void ewald::phys::calcSurfIntMat(NbrList &nlist, Mat lhs1=NULL, Mat lhs2=NULL)
{
    int nrow, ncol;

    // Argument check
    if (lhs1 == NULL && lhs2 == NULL) return;

    if (lhs1) MatGetSize(lhs1, &nrow, &ncol);
    if (lhs2) MatGetSize(lhs2, &nrow, &ncol);

    // Trivial case: zero-size matrix
    if (nrow == 0 || ncol == 0) {
      //        if (lhs1) {
      //      MatZeroEntries(lhs1);
      //      MatAssemblyBegin(lhs1, MAT_FINAL_ASSEMBLY);
      //      MatAssemblyEnd(lhs1, MAT_FINAL_ASSEMBLY);
      //  }

      //  if (lhs2) {
      //      MatZeroEntries(lhs2);
      //      MatAssemblyBegin(lhs2, MAT_FINAL_ASSEMBLY);
      //      MatAssemblyEnd(lhs2, MAT_FINAL_ASSEMBLY);
      //  }

        return;
    }

    // Create the matrix nonzero structures
    int *nnbrs = new int[nlist.verts.size()];
    int *nnz = new int[nrow];

    nlist.findNumNbrPoints(nnbrs);
    for (int i = 0; i < nlist.verts.size(); i++)
        nnz[3*i] = nnz[3*i+1] = nnz[3*i+2] = 3*nnbrs[i];

    // Allocate matrix storage
    if (lhs1) MatSeqAIJSetPreallocation(lhs1, 0, nnz);
    if (lhs2) MatSeqAIJSetPreallocation(lhs2, 0, nnz);

    delete []nnbrs;
    delete []nnz;

    // Compute matrix elements
    // Init
    if (lhs1) MatZeroEntries(lhs1);
    if (lhs2) MatZeroEntries(lhs2);

    for (int i = 0 ; i < nlist.verts.size(); i++) {
        Point *point = nlist.verts[i];

        for (int p = nlist.firstNbr[i]; p < nlist.firstNbr[i+1]; p++) {
            Tri *T = nlist.faces[p];

            double dist = nlist.dists[p];
            if (dist > rc) continue;

	    // Special case for double-layer potential
	    if ( lhs1 == NULL &&
	         ( point == T->vert[0] ||
		   point == T->vert[1] ||
		   point == T->vert[2] ) ) continue;

	    // Compute the interaction matrix
            double xtar[3];
            double lhs1_tmp[3][3][3], lhs2_tmp[3][3][3];
            double (*plhs1_tmp)[3][3], (*plhs2_tmp)[3][3];

            m_dcopy(3, point->x, xtar);
            to_CloseBy(T->xc, xtar);

            plhs1_tmp = lhs1 ? lhs1_tmp : NULL;
            plhs2_tmp = lhs2 ? lhs2_tmp : NULL;

            if (dist*dist > T->detJ)
                int_on_tri_regular(xtar, *T, plhs1_tmp, plhs2_tmp);
            else
                int_on_tri_singular(xtar, *T, plhs1_tmp, plhs2_tmp);

            // Map matrix elements to the global
            int irows[3], icols[9];

            for (int ii = 0; ii < 3; ++ii)
                irows[ii] = 3*i + ii;

            for (int l =  0; l < 3; ++l)
                for (int jj = 0; jj < 3; ++jj)
                    icols[3*l+jj] = 3*T->vert[l]->Gindx + jj;

            if (lhs1)
                MatSetValues(lhs1, 3, irows, 9, icols, &lhs1_tmp[0][0][0], ADD_VALUES);
            if (lhs2)
                MatSetValues(lhs2, 3, irows, 9, icols, &lhs2_tmp[0][0][0], ADD_VALUES);
        } // inbr
    } // i

    if (lhs1) {
        MatAssemblyBegin(lhs1, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(lhs1, MAT_FINAL_ASSEMBLY);
    }
    if (lhs2) {
        MatAssemblyBegin(lhs2, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(lhs2, MAT_FINAL_ASSEMBLY);
    }
}
