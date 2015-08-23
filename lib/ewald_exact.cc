#include "cxxheaders.h"
#include "mathfunc.h"
#include "ewald.h"

/* void ewald::calcSum
 * void ewald::calcPhysSum
 * void ewald::calcFourSum */


/* Calculate Ewald sum
 *     f(x) = sum_(y) f(y) G(y,x) + g(y) T(y,x) n(y)
 * Arguments:
 *  x -- target coordinates
 *  v -- target velocities (Ewald sum)
 *  y -- source coordinates
 *  f -- single-layer potential
 *  g -- double-layer potential
 *  nrml -- double-layer normal direction
 *  tol -- error toleratance */
void ewald::calcSum(const DoubArray2d &x, DoubArray2d &v,
	const DoubArray2d &y, 
	const DoubArray2d &f, const DoubArray2d &g, const DoubArray2d &nrml,
	double tol)
{
    const int Nx = x.size(0);
    v.resize(Nx,3);

    MArray<double,2> v1, v2;
    calcPhysSum(x, v1, y, f, g, nrml, tol);
    calcFourSum(x, v2, y, f, g, nrml, tol);

    v = v1;
    v += v2;
}


/* Analytically compute physical part of the Ewald sum */
void ewald::calcPhysSum(const DoubArray2d &x, DoubArray2d &v,
	const DoubArray2d &y, 
	const DoubArray2d &f, const DoubArray2d &g, const DoubArray2d &nrml,
	double tol)
{
    const int Nx = x.size(0);
    const int Ny = y.size(0);

    v.resize(Nx,3);
    v = 0.0;

    for (int ix = 0; ix < Nx; ix++)
    for (int iy = 0; iy < Ny; iy++) {
        double x0[3], y0[3], ftmp[3], gtmp[3], ntmp[3];
	FOR_D3 {
	    x0[d] = x(ix,d);
	    y0[d] = y(iy,d);

	    ftmp[d] = f(iy,d);
	    gtmp[d] = g(iy,d);
	    ntmp[d] = nrml(iy,d);
	}
	to_CloseBy(x0, y0);

	for (int n_max = 0; n_max < 10; n_max++) {
	    // Calculate the physical sum over "one-shell" at a time
	    double dv[3];
	    FOR_D3 dv[d] = 0.0;

	    for (int n0 = -n_max; n0 <= n_max; n0++)
	    for (int n1 = -n_max; n1 <= n_max; n1++)
	    for (int n2 = -n_max; n2 <= n_max; n2++) {

	        if ( n0 > -n_max && n0 < n_max &&
		     n1 > -n_max && n1 < n_max &&
		     n2 > -n_max && n2 < n_max ) continue;
		 
	        double xx[3];
		FOR_D3 xx[d] = y0[d] - x0[d];
		
		FOR_I3 {
		    xx[i] += n0*axis[0][i]
		           + n1*axis[1][i]
			   + n2*axis[2][i];
		}

		double rr = normalizeVec3D(xx);
		double S1, S2, D1;
		phys::coeff_exact(rr, &S1, &S2, &D1);

		double xxf = m_ddot(3, xx, ftmp);
		double xxg = m_ddot(3, xx, gtmp);
		double xxn = m_ddot(3, xx, ntmp);

		m_daxpy(3, S1*xxf, xx, dv);
		m_daxpy(3, S2, ftmp, dv);
		m_daxpy(3, D1*xxg*xxn, xx, dv);
	    }

	    // Add the contribution from the shell
	    m_dadd(3, dv, &v(ix,0));

	    double dv2 = normVec3D(dv);
	    double v2 = normVec3D(&v(ix,0));
	    if (dv2 < (v2 + 1.E-6) *tol) break;;
	} // n_max
    } // ix, iy
}


/* Compute Fourier part of the Ewald sum exactly */
void ewald::calcFourSum(const DoubArray2d &x, DoubArray2d &v,
	const DoubArray2d &y, 
	const DoubArray2d &f, const DoubArray2d &g, const DoubArray2d &nrml,
	double tol)
{
    const int Nx = x.size(0);
    const int Ny = y.size(0);

    const dcomplex iota(0.0, 1.0);
    const double iL[3] = { 1/L[0], 1/L[1], 1/L[2] };
    const double vol = L[0]*L[1]*L[2];

    // Init
    v.resize(Nx,3);
    v = 0.0;

    for (int n_max = 1; n_max < 10000; n_max++) {

	// Inverse transform back to target points
	MArray<double,2> dv(Nx,3);
	dv = 0.0;

	for (int n0 = -n_max; n0 <= n_max; n0++)
	for (int n1 = -n_max; n1 <= n_max; n1++)
	for (int n2 = -n_max; n2 <= n_max; n2++) {
	    if ( n0 > -n_max && n0 < n_max &&
		 n1 > -n_max && n1 < n_max &&
		 n2 > -n_max && n2 < n_max ) continue;

	    double q[3], qt[3], q2, q2t, expq2t, phi0, phi1;
	    
	    FOR_I3 {
	        q[i] = n0*rcp_axis[0][i]
		     + n1*rcp_axis[1][i]
		     + n2*rcp_axis[2][i];
	    }

	    FOR_D3 qt[d] = sqrt(M_PI*alpha)*q[d];
            q2 = m_ddot(3, q, q);
            q2t = M_PI*alpha*q2;
            expq2t = exp(-q2t);

            phi0 = expq2t/q2t;
            phi1 = (expq2t + phi0)/q2t;

            double S1 = 2*alpha/vol*phi1;
	    dcomplex D1 = 4*M_PI*alpha/vol*phi0*iota;
	    dcomplex D2 = -8*(M_PI*M_PI)*(alpha*alpha)/vol*phi1*iota;
	    
	    // Fourier transform the source term
	    dcomplex coeff[3];
	    FOR_D3 coeff[d] = 0.0;

	    for (int iy = 0; iy < Ny; iy++) {
		double ytmp[3], ftmp[3], gtmp[3], ntmp[3];

		FOR_D3 {
		    ytmp[d] = y(iy,d);
		    ftmp[d] = f(iy,d);
		    gtmp[d] = g(iy,d);
		    ntmp[d] = nrml(iy,d);
	        }

		double qy = m_ddot(3, q, ytmp);
		dcomplex expiqy = exp(2*M_PI*qy*iota);

		// Single-layer
		{
		    double qtf = m_ddot(3, qt, ftmp);

		    FOR_D3 coeff[d] += 
		    	S1*(q2t*ftmp[d] - qt[d]*qtf)*expiqy;
	        }

		// Double-layer
		{
                    double gn = m_ddot(3, gtmp, ntmp);
		    double gq = m_ddot(3, gtmp, q);
		    double qn = m_ddot(3, q, ntmp);

		    FOR_D3 {
		        coeff[d] += 
				( D1*(gn*q[d] + gq*ntmp[d] + qn*gtmp[d]) 
				+ D2*gq*qn*q[d] )*expiqy;
                    }
                }
	    } // iy

	    for (int ix = 0; ix < Nx; ix++) {
	        double qx = m_ddot(3, q, &x(ix,0));
		dcomplex expiqx = exp(-2*M_PI*qx*iota);
		FOR_D3 dv(ix,d) += real(expiqx*coeff[d]);
	    }
	} // n0, n1, n2

	v += dv;

	double v2 = 0.0, dv2 =  0.0;
	for (int ix = 0; ix < Nx; ix++) {
	    FOR_D3 {
		v2 += v(ix,d)*v(ix,d);
		dv2 += dv(ix,d)*dv(ix,d);
	    }
	}
	v2 = sqrt(v2);
	dv2 = sqrt(dv2);

	if (dv2 < (v2 + 1.E-6)*tol) break;
    } // n_max
}
