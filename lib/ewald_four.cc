#include "cxxheaders.h"
#include "mathfunc.h"
#include "ewald.h"
#include "tri.h"
#include "rfftw_mpi.h"

/* bool ewald::four::point_active 
 * bool ewald::four::tri_active 
 * void ewald::four::getBlkNum
 * void ewald::four::setSources
 * void ewald::four::setTargets
 * void ewald::four::initPme
 * void ewald::four::clear_source
 * void ewald::four::add_source
 * void ewald::four::add_interp_vel
 * void ewald::four::transform */


bool ewald::four::point_active(const double *x)
{
    double xtmp = (ewald::bctype == RECT) ?  x[0] :
		L[0]*( rcp_axis[0][0]*x[0] + rcp_axis[0][2]*x[2] );
    			
    return point_in_seg(xtmp, xminBuf, xmaxBuf, ewald::L[0]);
}


/* Whether a triangle can have non-zero contribution to local
 * physical Ewald sum
 * Arguments:
 *   x[i] -- the i-th point */
bool ewald::four::tri_active(const double (*x)[3])
{
    double xminTri = FLT_MAX;
    double xmaxTri = -FLT_MAX;

    for (int l = 0; l < 3; l++) {
        double xtmp = (ewald::bctype == RECT) ?  x[l][0] :
		L[0]*( rcp_axis[0][0]*x[l][0] + rcp_axis[0][2]*x[l][2] );

        xminTri = min(xminTri, xtmp);
	xmaxTri = max(xmaxTri, xtmp);
    }

    return seg_overlap(xminBuf, xmaxBuf, xminTri, xmaxTri, ewald::L[0]);
}


/* Find which block the a point lies in
 * Arguments:
 *  x -- point coordinate
 *  xb -- block number
 * Note:
 *   iBgn[0] <= xb[0] < iEnd[0] */
void ewald::four::getBlkNum(const double *x, double *xb)
{
    phys_to_lattice_coord(x, xb);
    FOR_D3 xb[d] *= Nb[d];
}


/* Find sources
 * Arguments:
 *   allSrcs -- all available sources
 *   srcs -- active sources */
void ewald::four::setSources(const vector<Tri*> &tris, vector<Tri*> &srcs)
{
    srcs.clear();

    for (int i = 0; i < tris.size(); i++) {
        Tri *T = tris[i];
	double x[3][3];

	FOR_J3
	FOR_D3 {
	    x[j][d] = T->vert[j]->x[d];
	}

	if (tri_active(x)) srcs.push_back(T);
    }
}


/* Find targets
 * Arguments:
 *   points -- all possible targets
 *   tgts -- active tagets  */
void ewald::four::setTargets(const vector<Point*> &points, vector<Point*> &tgts)
{
    tgts.clear();

    for (int i = 0; i < points.size(); i++) {
        Point *vert = points[i];
	if (point_active(vert->x)) tgts.push_back(vert);
    }
}


/* Initialie the PME solver
 * Note:
 *   1. create FFTW plans
 *   2. allocate memories 
 *   3. pre-compute transform coefficients */ 
void ewald::four::initPme()
{
    assert(four::active);

    // Domain decomposition
    xmin = comm_rank*L[0]/comm_size;
    xmax = (comm_rank + 1)*L[0]/comm_size;

    xminBuf = xmin;
    xmaxBuf = xmax + PB*L[0]/Nb[0];

    iBgn[0] = comm_rank*Nb[0]/comm_size;
    iEnd[0] = iBgn[0] + Nb[0]/comm_size;
    iBgn[1] = 0;    iEnd[1] = Nb[1];
    iBgn[2] = 0;    iEnd[2] = Nb[2];
    
    qBgn[0] = 0;    qEnd[0] = Nb[0];
    qBgn[1] = comm_rank*Nb[1]/comm_size;
    qEnd[1] = qBgn[1] + Nb[1]/comm_size;
    qBgn[2] = 0;    qEnd[2] = Nb[2]/2+1;

    // FFTW2 plans
    _fplan = rfftw3d_mpi_create_plan(comm, Nb[0], Nb[1], Nb[2],
		  FFTW_REAL_TO_COMPLEX, FFTW_ESTIMATE);

    _bplan = rfftw3d_mpi_create_plan(comm, Nb[0], Nb[1], Nb[2],
		  FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE);
    
    { // check domain decomposition is the same as required by FFTW2
	int local_nx, local_x_start, local_ny_after_transpose,
	      local_y_start_after_transpose, total_local_size;


	rfftwnd_mpi_local_sizes(_fplan, &local_nx, &local_x_start, 
		      &local_ny_after_transpose,
		      &local_y_start_after_transpose,
		      &total_local_size);

	assert(local_nx == Nb[0]/comm_size);
	assert(local_x_start == iBgn[0]);
	assert(local_ny_after_transpose == Nb[1]/comm_size);
	assert(local_y_start_after_transpose == qBgn[1]);
    }


    {  
	IndexRange r0(iBgn[0], iEnd[0]);
	IndexRange r1(iBgn[1], iEnd[1]);
	IndexRange r2(iBgn[2], iEnd[2]);

        _f.resize(r0,r1,r2,3);
	_g.resize(r0,r1,r2,9);
	_v.resize(r0,r1,r2,3);
    }

    { // compute b-factor
        const dcomplex iota(0, 1);
	double MP[PB];

	// allocate memory for bb
	bb0.resize(Nb[0]);
	bb1.resize(Nb[1]);
	bb2.resize(Nb[2]);

	int imin;
	computeBSpline(double(PB), PB, imin, MP);
	if (imin /= 1) {
	    computeBSpline(double(PB)+DBL_EPSILON, PB, imin, MP);
	}

	FOR_D3 {
	    for (int k = qBgn[d]; k < qEnd[d]; ++k) {
		double c;
	        dcomplex b = 0;

	        for (int m = 0; m < PB-1; m++) {
		    c = 2*M_PI*k*m/double(Nb[d]);
	            b += MP[m]*exp(c*iota);
		}

		c = 2*M_PI*(PB-1.0)/double(Nb[d]);
		b = exp(c*iota)/b;
		b = norm(b);

		switch (d) {
		    case 0:
			bb0(k) = b.real();
			break;
		    case 1:
			bb1(k) = b.real();
			break;
		    case 2:
			bb2(k) = b.real();
			break;
		}
	    } // k
        } // d
    }
}


/* Zero out sources and reset single and double layer density flags
 * i.e. _f_flag and _g_flag */
void ewald::four::clear_source()
{
    _f = 0.0;
    _g = 0.0;
    _v = 0.0;

    _f_flag = false;
    _g_flag = false;
}


/*  Distribute point sources
 *  Arguments:
 *    n -- number of singularities
 *    x -- coordinates of singularities
 *    cf -- coefficients before the single-layer density
 *    f -- single-layer singularities
 *    cg -- coefficients before the double-layer density
 *    g -- double layer singularity
 *    normal -- the surface normal */
void ewald::four::add_source(int n, const double *x, 
		double cf, const double *f, 
		double cg, const double *g, const double *normal)
{
    if (n <= 0) return;

    // Update flags
    bool _has_f = fabs(cf) > 1.E-10;
    bool _has_g = fabs(cg) > 1.E-10;

    _f_flag = _f_flag || _has_f;
    _g_flag = _g_flag || _has_g;
    
    for (int i = 0; i < n; ++i) {
        double xtmp[3], ftmp[3], gtmp[3], ntmp[3];
	double xb[3];
	int ix0, iy0, iz0, ix, iy, iz;
	double wx[PB], wy[PB], wz[PB], wxyz;

	m_dcopy(3, &x[3*i], xtmp);
	if (_has_f) 
	    m_dcopy(3, &f[3*i], ftmp);
	if (_has_g) {
	    m_dcopy(3, &g[3*i], gtmp);
	    m_dcopy(3, &normal[3*i], ntmp);
        }

	if (!ewald::four::point_active(xtmp)) continue;

	getBlkNum(xtmp, xb);
	computeBSpline(xb[0], PB, ix0, wx);
	computeBSpline(xb[1], PB, iy0, wy);
	computeBSpline(xb[2], PB, iz0, wz);

	for (int dix = 0; dix < PB; dix++) {
	    ix = modulo(ix0 + dix, Nb[0]);
	    if (ix < iBgn[0] || ix >= iEnd[0]) continue;

	    for (int diy = 0; diy < PB; diy++)
	    for (int diz = 0; diz < PB; diz++)  {
		iy = modulo(iy0 + diy, Nb[1]); 
		iz = modulo(iz0 + diz, Nb[2]);

		wxyz = wx[dix] * wy[diy] * wz[diz];

		if (_has_f) {
		    for (int ii = 0; ii < 3; ++ii)
		      _f(ix,iy,iz,ii) += wxyz * cf * ftmp[ii];
		}

		if (_has_g) {
		    for (int ii = 0; ii < 3; ++ii)
		    for (int jj = 0; jj < 3; ++jj)
		      _g(ix,iy,iz,3*ii+jj) += wxyz * cg * gtmp[ii] * ntmp[jj];
		}
	    } // diy, diz
	} // dix
    } // i
}


/* Add source from a source list
 * Arguments:
 *   cf, cg -- coefficients of the single- and double-layer potentials
 *   slist -- the source list */
void ewald::four::add_source(double cf, double cg, const vector<Tri*> &slist)
{
    // Set flags
    bool _has_f = fabs(cf) > 1.E-10;
    bool _has_g = fabs(cg) > 1.E-10;
    
    // Use 3-point quadrature rule
    Quad2D &Q = quadrature::select_rule_2d("TRI_3");
    int nq = Q.n();
    MArray<double,2> xq(nq,3), fq(nq,3), gq(nq,3), normalq(nq,3);

    for (int i = 0; i < slist.size(); ++i) {
        Tri *face = slist[i];
        double xtri[3][3], ftri[3][3], gtri[3][3];

	for (int l = 0; l < 3; ++l) {
	    m_dcopy(3, face->vert[l]->x, xtri[l]);

	    if (_has_f) m_dcopy(3, face->vert[l]->f, ftri[l]);
            if (_has_g) m_dcopy(3, face->vert[l]->g, gtri[l]);
	}

	for (int iq = 0; iq < Q.n(); ++iq) {
	    double s = Q.x(iq), t = Q.y(iq);
	    double dA = face->detJ*Q.w(iq);
	    double w0, w1, w2;		// weights
	    w0 = 1.0 - s - t;
	    w1 = s;
	    w2 = t;

	    for (int ii = 0; ii < 3; ++ii) {
		xq(iq,ii) = w0*xtri[0][ii] + w1*xtri[1][ii] + w2*xtri[2][ii];

		if (_has_f) {
		    fq(iq,ii) = w0*ftri[0][ii] + w1*ftri[1][ii] + w2*ftri[2][ii];
		    fq(iq,ii) *= dA;
		} else {
		    fq(iq,ii) = 0.0;
                }

		if (_has_g) {
		    gq(iq,ii) = w0*gtri[0][ii] + w1*gtri[1][ii] + w2*gtri[2][ii];
		    gq(iq,ii) *= dA;
		    normalq(iq,ii) = face->normal[ii];
	        } else {
		    gq(iq,ii) = 0.0;
		    normalq(iq,ii) = 0.0;
		}
	    } // i
	} // iq

	add_source(nq, xq.data(), cf, fq.data(), cg, gq.data(), normalq.data());
    }  // i
}


/* Add the interpolated velocity from the mesh
 * Arguments:
 *   Nx -- number of target points
 *   x -- target point coordinates
 *   v -- target point velocities */
void ewald::four::add_interp_vel(int Nx, const double *x, double *v)
{
    for (int i = 0; i < Nx; ++i) {
	double xtmp[3], xb[3];
	double wx[PB], wy[PB], wz[PB], wxyz;
	int ix0, iy0, iz0, ix, iy, iz;

	m_dcopy(3, &x[3*i], xtmp);

	getBlkNum(xtmp, xb);
	computeBSpline(xb[0], PB, ix0, wx);
	computeBSpline(xb[1], PB, iy0, wy);
	computeBSpline(xb[2], PB, iz0, wz);

	for (int dix = 0; dix < PB; dix++) {
	    ix = modulo(ix0 + dix, Nb[0]);
	    if (ix < iBgn[0] || ix >= iEnd[0]) continue;

	    for (int diy = 0; diy < PB; diy++) 
	    for (int diz = 0; diz < PB; diz++)  {
		iy = modulo(iy0 + diy, Nb[1]);
		iz = modulo(iz0 + diz, Nb[2]);

		wxyz = wx[dix] * wy[diy] * wz[diz];

		for (int d = 0; d < 3; ++d)
		    v[3*i+d] += wxyz*_v(ix,iy,iz,d);
	    } // diy, diz
	} // dix
    } // i
}


// Overloaded add_interp_vel
void ewald::four::add_interp_vel(const vector<Point*> &tgts, double *v)
{
    int n = tgts.size();
    double *x = new double[3*n];

    for (int i = 0; i < n; ++i) {
        m_dcopy(3, tgts[i]->x, x+3*i);
    }

    add_interp_vel(n, x, v);

    // delete working arrays
    delete [] x;
}


// Do the transform on the Cartesian mesh
void ewald::four::transform()
{
    int n0, n1, n2, nvar, lwork;
    int stride_x[3], stride_q[3];
    double *data, *work;
    fftw_complex *cdata;

    // Synchronize _f_flag and _g_flag
    int _my_f_flag = _f_flag? 1 : 0; 
    int _my_g_flag = _g_flag? 1 : 0;
    int _all_f_flag, _all_g_flag;

    MPI_Allreduce(&_my_f_flag, &_all_f_flag, 1, MPI_INT, MPI_MAX, comm);
    MPI_Allreduce(&_my_g_flag, &_all_g_flag, 1, MPI_INT, MPI_MAX, comm);

    _f_flag = _all_f_flag;
    _g_flag = _all_g_flag;


    // Trivial case
    if ((!_f_flag) && (!_g_flag)) {
        _v = 0.0;
	return;
    }

    // Allocate memory
    // Assuming nx and ny are divisible by comm_size
    nvar = 0;
    if (_f_flag) nvar += 3;
    if (_g_flag) nvar += 9;

    n0 = Nb[0];
    n1 = Nb[1];
    n2 = Nb[2];
    lwork = n0/comm_size * n1 * 2*(n2/2 + 1);
    lwork *= nvar;

    data = new double[lwork];
    cdata = (fftw_complex*)data;
    work = new double[lwork];

    // calculate memory stride
    stride_x[2] = 1;
    stride_x[1] = stride_x[2] * 2*(n2/2+1);
    stride_x[0] = stride_x[1] * n1;

    // note: the array is transposed by FFTW
    // for Foruier coefficient, the order of index is (y, x, z) for
    // index changes from the fastest to slowest
    stride_q[2] = 1;
    stride_q[0] = stride_q[2] * (n2/2+1);
    stride_q[1] = stride_q[0] * n0;

    // copy data
    for (int i0 = iBgn[0]; i0 < iEnd[0]; ++i0)
    for (int i1 = iBgn[1]; i1 < iEnd[1]; ++i1)
    for (int i2 = iBgn[2]; i2 < iEnd[2]; ++i2) {
	int offset = i2 + i1*stride_x[1] + (i0 - iBgn[0])*stride_x[0];
	offset *= nvar;

	if (_f_flag) {
	    m_dcopy(3, &_f(i0,i1,i2,0), data+offset);
	    offset += 3;
	}

	if (_g_flag) {
	    m_dcopy(9, &_g(i0,i1,i2,0), data+offset);
	    offset += 9;
	}
    } // i0 i1 i2

    // forward FFT
    rfftwnd_mpi(_fplan, nvar, data, work, FFTW_TRANSPOSED_ORDER);

    // Transform the Fourier coefficients
    for (int i1 = qBgn[1]; i1 < qEnd[1]; i1++)
    for (int i0 = qBgn[0]; i0 < qEnd[0]; i0++)
    for (int i2 = qBgn[2]; i2 < qEnd[2]; i2++) {
        // q -- wave number
	// qlat -- wave number in the reciprocal lattice space
        int qlat[3];		
	double q[3];

        qlat[0] = i0 < n0/2 ? i0 : (i0 - n0);
        qlat[1] = i1 < n1/2 ? i1 : (i1 - n1);
        qlat[2] = i2;

	q[0] = qlat[0]*rcp_axis[0][0] + qlat[2]*rcp_axis[2][0];
	q[1] = qlat[1]*rcp_axis[1][1];
	q[2] = qlat[0]*rcp_axis[0][2] + qlat[2]*rcp_axis[2][2];

	q[0] = -q[0];  // Set notes for why q should be set to -q
	q[1] = -q[1];
	q[2] = -q[2];

	// copy Fourier component
        dcomplex fh[3], gh[3][3], vh[3];

        long offset = i2 + i0*stride_q[0] + (i1 - qBgn[1])*stride_q[1];
        offset *= nvar;
	fftw_complex *p = cdata + offset;

        if (_f_flag) {
	    for (int l = 0; l < 3; ++l) {
	        fh[l] = dcomplex(p->re, p->im);
	        ++p;
	    }
        }

	if (_g_flag) {
	    for (int l = 0; l < 3; ++l)
	    for (int m = 0; m < 3; ++m) {
	        gh[l][m] = dcomplex(p->re, p->im);
		++p;
	    }
	}

	// Transform Fourier coefficients
	vh[0] = vh[1] = vh[2] = 0.0;

	if (i0 != 0 || i1 != 0 || i2 != 0) {
	    double qt[3], q2, q2t, expq2t, phi0, phi1;
	    double A;

	    A = sqrt(M_PI*alpha);
	    for (int d = 0; d < 3; ++d) qt[d] = A * q[d];

	    q2 = q[0]*q[0] + q[1]*q[1] + q[2]*q[2];
	    q2t = (M_PI*alpha)*q2;
	    expq2t = exp(-q2t);

	    phi0 = expq2t/q2t;
	    phi1 = (expq2t + phi0)/q2t;

	    if (_f_flag) {
	        dcomplex qtf = qt[0]*fh[0] + qt[1]*fh[1] + qt[2]*fh[2];
		double A = 2*alpha*iVol*phi1;
		for (int l = 0; l < 3; ++l)
		    vh[l] += A*(q2t*fh[l] - qt[l]*qtf);
            }

	    if (_g_flag) {
	        dcomplex trg = gh[0][0] + gh[1][1] + gh[2][2];
		dcomplex qg[3], gq[3], qgq;

		for (int l = 0; l < 3; ++l) {
		    qg[l] = q[0]*gh[0][l] + q[1]*gh[1][l] + q[2]*gh[2][l];
		    gq[l] = gh[l][0]*q[0] + gh[l][1]*q[1] + gh[l][2]*q[2];
                }
		qgq = q[0]*gq[0] + q[1]*gq[1] + q[2]*gq[2];

		const dcomplex iota(0.0, 1.0);
		double A1 = 4*M_PI*alpha*iVol*phi0;
		double A2 = 8*(M_PI*M_PI)*(alpha*alpha)*iVol*phi1;
		
		for (int l = 0; l < 3; ++l) {
		    vh[l] += iota*A1*(q[l]*trg + qg[l] + gq[l]) - iota*A2*qgq*q[l];
		}
	    }

	    // multiply vh by b-factor
	    double bb = bb0(i0) * bb1(i1) * bb2(i2);
	    for (int l = 0; l < 3; l++) vh[l] *= bb;
        }


	// copy vh back to data
	{
            long offset = i2 + i0*stride_q[0] + (i1 - qBgn[1])*stride_q[1];
	    offset *= 3;
	    fftw_complex *p = cdata + offset;

	    for (int d = 0; d < 3; ++d) {
	        p->re = vh[d].real();
	        p->im = vh[d].imag();
		p++;
	    }
        }
    } // i1,2,3

    // backward FFTW
    rfftwnd_mpi(_bplan, 3, data, work, FFTW_TRANSPOSED_ORDER);

    // copy data back
    for (int i0 = iBgn[0]; i0 < iEnd[0]; ++i0)
    for (int i1 = iBgn[1]; i1 < iEnd[1]; ++i1)
    for (int i2 = iBgn[2]; i2 < iEnd[2]; ++i2) {
	int offset = i2 + i1*stride_x[1] + (i0 - iBgn[0])*stride_x[0];
	offset *= 3;
	m_dcopy(3, data+offset, &_v(i0,i1,i2,0));
    }

    // deallocate working arrays
    delete []data;
    delete []work;
}
