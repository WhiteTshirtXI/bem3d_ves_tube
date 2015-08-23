#include "cxxheaders.h"
#include "subdiv.h"
#include "mathfunc.h"
#include "mesh.h"
#include "quadrature.h"
#include "marray.h"

/* void subdiv::divideMesh
 * void subdiv::calcFunc_666 
 * void subdiv::calcFunc_Recursive
 * void subdiv::calcFunc
 * void subdiv::buildFuncDataBase
 * void subdiv::buildControlPoints */


/* Apply Loop subdivision once 
 * Arguments:
 *  x, f2v -- the mesh coordinates and the connectivity
 *  xnew, f2vnew -- the coordinates and connectivity after subidvision */
void subdiv::divideMesh(const MArray<double,2> &x, const MArray<int,2> &f2v,
                MArray<double,2> &x_new, MArray<int,2> &f2v_new)
{
    int nvert = x.size(0);
    int nface = f2v.size(0);

    // Build edge ist
    int nedge_max = 3*nface;
    int (*e2v)[2] = new int[nedge_max][2];
    int (*e2f)[2] = new int[nedge_max][2];
    int (*f2e)[3] = new int[nface][3];
    int (*f2f)[3] = new int[nface][3];
    int nedge;

    {
        int (*f2v_tmp)[3] = new int[nface][3];

	for (int iface = 0; iface < nface; iface++)
	for (int l = 0; l < 3; l++) {
	    f2v_tmp[iface][l] = f2v(iface,l);
	}

	Mesh_edgeList(nvert, nface, f2v_tmp, nedge, e2v, e2f);
	Mesh_faceNbrList(nvert, nface, f2v_tmp, nedge, e2v, e2f, f2e, f2f);
	delete [] f2v_tmp;
    }

    // Calculate the new mesh coordinates
    MArray<double,2> xnode(nvert,3), xedge(nedge,3);
    xnode = 0.0;
    xedge = 0.0;

    // Edge points
    for (int iedge = 0; iedge < nedge; iedge++) {
        int i0 = e2v[iedge][0];
	int i1 = e2v[iedge][1];

	m_daxpy(3, 0.375, &x(i0,0), &xedge(iedge,0));
	m_daxpy(3, 0.375, &x(i1,0), &xedge(iedge,0));

	for (int l = 0; l < 2; l++) {
	    int iface = e2f[iedge][l];
	    for (int m = 0; m < 3; m++) {
	        int i = f2v(iface,m);
		if (i != i0 && i != i1) {
		    m_daxpy(3, 0.125, &x(i,0), &xedge(iedge,0));
		    break;
		}
	    }
	}
    }

    // Nodal points 
    {
	// Find the number of points
        MArray<int,1> nnbr(nvert);
	nnbr = 0;

	for (int iedge = 0; iedge < nedge; iedge++) {
	    int i0 = e2v[iedge][0];
	    int i1 = e2v[iedge][1];

	    nnbr(i0)++;
	    nnbr(i1)++;
        }

	MArray<double,1> w(nvert);
	for (int ivert = 0; ivert < nvert; ivert++) {
	    double tmp = 0.375 + 0.25*cos(2*M_PI/nnbr(ivert));
	    w(ivert) = 1.0/nnbr(ivert)*(0.675 - tmp*tmp);

	    m_daxpy(3, 1.0-nnbr(ivert)*w(ivert), &x(ivert,0), &xnode(ivert,0));
	}

	for (int iedge = 0; iedge < nedge; iedge++) {
	    int i0 = e2v[iedge][0];
	    int i1 = e2v[iedge][1];

	    m_daxpy(3, w(i0), &x(i1,0), &xnode(i0,0));
	    m_daxpy(3, w(i1), &x(i0,0), &xnode(i1,0));
        }
    }


    // Now fill the entries
    x_new.resize(nvert+nedge, 3);
    m_dcopy(3*nvert, xnode.data(), x_new.data());
    m_dcopy(3*nedge, xedge.data(), x_new.data() + 3*nvert );

    f2v_new.resize(4*nface, 3);
    for (int cnt = 0, iface = 0; iface < nface; iface++) {
        int i0 = f2v(iface,0);
        int i1 = f2v(iface,1);
        int i2 = f2v(iface,2);

	int j0 = f2e[iface][0];
	int j1 = f2e[iface][1];
	int j2 = f2e[iface][2];

	f2v_new(cnt,0) = i0;
	f2v_new(cnt,1) = j2 + nvert;
	f2v_new(cnt,2) = j1 + nvert;
	cnt++;

	f2v_new(cnt,0) = i1;
	f2v_new(cnt,1) = j0 + nvert;
	f2v_new(cnt,2) = j2 + nvert;
	cnt++;

	f2v_new(cnt,0) = i2;
	f2v_new(cnt,1) = j1 + nvert;
	f2v_new(cnt,2) = j0 + nvert;
	cnt++;

	f2v_new(cnt,0) = j0 + nvert;
	f2v_new(cnt,1) = j1 + nvert;
	f2v_new(cnt,2) = j2 + nvert;
	cnt++;
    }

    // Delete temp arrays
    delete [] e2v;
    delete [] e2f;
    delete [] f2e;
    delete [] f2f;
}


/* Compute regular Loop subdivision shape function
 * Arguments:
 *  v, w -- the barycentric coordinates
 *  f -- function value
 *  df -- df/dv, df/dw
 *  ddf -- ddf/dvdv ddf/dwdw ddf/dvdw */
void subdiv::calcFunc_666(double v, double w, 
	MArray<double,1> &f, MArray<double,2> &df, MArray<double,3> &ddf)
{
    const double phi[12][15] = {
	{ 6, 0, 0, -12, -12, -12, 8, 12, 12, 8, -1, -2, 0, -2, -1 },
	{ 1, 4, 2, 6, 6, 0, -4, -6, -12, -4, -1, -2, 0, 4, 2},
	{ 1, 2, 4, 0, 6, 6, -4, -12, -6, -4, 2, 4, 0, -2, -1},
	{ 0, 0, 0, 0, 0, 0, 2, 6, 6, 2, -1, -2, 0, -2, -1},
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1},
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, -1},
	{ 1, -2, 2, 0, -6, 0, 2, 6, 0, -4, -1, -2, 0, 4, 2},
	{ 1, -4, -2, 6, 6, 0, -4, -6, 0, 2, 1, 2, 0, -2, -1},
	{ 1, -2, -4, 0, 6, 6, 2, 0, -6, -4, -1, -2, 0, 2, 1},
	{ 1, 2, -2, 0, -6, 0, -4, 0, 6, 2, 2, 4, 0, -2, -1},
	{ 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -1, -2, 0, 0, 0},
	{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0}
    };

    f.resize(12);
    df.resize(12,2);
    ddf.resize(12,2,2);

    // f
    { 
        f = 0.0;

        double c[] = { 1, 
			v, w, 
			v*v, v*w, w*w, 
			v*v*v, v*v*w, v*w*w, w*w*w,
			v*v*v*v, v*v*v*w, v*v*w*w, v*w*w*w, w*w*w*w };

        for (int i = 0; i < 12; i++) {
	    for (int j = 0; j < 15; j++) {
	        f(i) += 1.0/12*phi[i][j]*c[j];
	    }
	}
    }

    { // df
	df = 0.0;

        double dv[] = { 0, 
			1, 0, 
			2*v, w, 0,
			3*v*v, 2*v*w, w*w, 0,
			4*v*v*v, 3*v*v*w, 2*v*w*w, w*w*w, 0};
        double dw[] = { 0, 
			0, 1, 
			0, v, 2*w,
			0, v*v, 2*v*w, 3*w*w,
			0, v*v*v, 2*v*v*w, 3*v*w*w, 4*w*w*w};

        for (int i = 0; i < 12; i++) {
	    for (int j = 0; j < 15; j++) {
		df(i,0) += 1.0/12*phi[i][j]*dv[j];
		df(i,1) += 1.0/12*phi[i][j]*dw[j];
	    }
	}
    }

    { // ddf
	ddf = 0.0;

        double dvv[] = { 0, 
			0, 0,
			2, 0, 0,
			6*v, 2*w, 0, 0,
			12*v*v, 6*v*w, 2*w*w, 0, 0};
        double dww[] = { 0, 
			0, 0,
			0, 0, 2,
			0, 0, 2*v, 6*w,
			0, 0, 2*v*v, 6*v*w, 12*w*w };
        double dvw[] = { 0, 
			0, 0,
			0, 1, 0,
			0, 2*v, 2*w, 0,
			0, 3*v*v, 6*v*w, 3*w*w, 0 };

        for (int i = 0; i < 12; i++) {
	    for (int j = 0; j < 15; j++) {
	        ddf(i,0,0) += 1.0/12*phi[i][j]*dvv[j];
		ddf(i,1,1) += 1.0/12*phi[i][j]*dww[j];
		ddf(i,0,1) += 1.0/12*phi[i][j]*dvw[j];
	    }

            ddf(i,1,0) = ddf(i,0,1);	// symmetry
	}
    }
}


/* Calculate subdivision function
 * Arguments:
 *  N0,N1,N2 -- valance of the three vertices
 *  v,w -- barycentric coordinate
 *  f, df, ddf -- 
 *  level -- level of subdivision 
 *  min_level -- minimum level of subdivision, for debug only */
void subdiv::calcFunc_Recursive(int N0, int N1, int N2, double v, double w, 
	MArray<double,1> &f, MArray<double,2> &df, MArray<double,3> &ddf,
	int level, int min_level)
{
    assert(N0 > 3 && N1 > 3 && N2 > 3);

    const int MAX_LEVEL = 25;

    if (N0 == 6 && N1 == 6 && N2 == 6 && level >= min_level) {
	if (min_level >= 0) printf("Subdivision stops at level %2d\n", level);

        calcFunc_666(v, w, f, df, ddf);
	return;
    }

    // Move the point to the center after too many subdivisions
    if (level >= MAX_LEVEL) {
        v = w = 1.0/3;
    }

    // Subdivde the mesh
    int N = N0 + N1 + N2 - 6;
    double u = 1 - v - w;
    int case_number;
    int N0_fn, N1_fn, N2_fn, N_fn;
    double v_fn, w_fn;
    MArray<double,1> f_fn;
    MArray<double,2> df_fn;
    MArray<double,3> ddf_fn;

    if (u > 0.5) {
	case_number = 0;
	N0_fn = N0;
	N1_fn = 6;
	N2_fn = 6;
	v_fn = 2*v;
	w_fn = 2*w;
    }
    else if (v > 0.5) {
	case_number = 1;
	N0_fn = 6;
	N1_fn = N1;
	N2_fn = 6;
	v_fn = 2*v - 1;
	w_fn = 2*w;
    }
    else if (w > 0.5) {
	case_number = 2;
	N0_fn = 6;
	N1_fn = 6;
	N2_fn = N2;
	v_fn = 2*v;
	w_fn = 2*w - 1;
    }
    else {
	case_number = 3;
	N0_fn = 6;
	N1_fn = 6;
	N2_fn = 6;
	v_fn = 1 - 2*v;
	w_fn = 1 - 2*w;
    }

    N_fn = N0_fn + N1_fn + N2_fn - 6;
    calcFunc_Recursive(N0_fn, N1_fn, N2_fn, v_fn, w_fn, f_fn, df_fn, ddf_fn, level+1, min_level);

    // Now map the functions back to the original mesh
    // First find the neighbors of the three vertices 
    vector<int> nbrs[3];
    {
	int p = 3;
	nbrs[2].push_back(1);
	nbrs[2].push_back(p);
	for (int i = 1; i < N2-2; i++) nbrs[2].push_back(++p);
	nbrs[2].push_back(0);

	nbrs[0].push_back(2);
	nbrs[0].push_back(p);
	for (int i = 1; i < N0-2; i++) nbrs[0].push_back(++p);
	nbrs[0].push_back(1);

	nbrs[1].push_back(0);
	nbrs[1].push_back(p);
	for (int i = 1; i < N1-3; i++) nbrs[1].push_back(++p);
	nbrs[1].push_back(3);
	nbrs[1].push_back(2);
    }


    // Label each point of the finer mesh with two integers.
    // -- If the point is the i-th mesh point on the coarser mesh, we label it as (i,i).
    // -- If the point is the mid-point of an edge of the coarser mesh, we label it
    //    as (i,j), where i, j are the two end points of that edge.
    MArray<int,2> label(N_fn,2);

    if (case_number == 0) {
        int p = 0;
	label(p,0) = 0;    label(p++,1) = 0;
	label(p,0) = 0;    label(p++,1) = 1;
	label(p,0) = 0;    label(p++,1) = 2;

	label(p,0) = 1;    label(p++,1) = 2;
	label(p,0) = 2;    label(p++,1) = 2;
	label(p,0) = 2;    label(p++,1) = nbrs[2][N2-2];
	for (int i = 1; i < N0-1; i++) {
	    label(p,0) = 0;    label(p++,1) = nbrs[0][i];
	}
	label(p,0) = 1;    label(p++,1) = nbrs[1][1];
	label(p,0) = 1;    label(p++,1) = 1;

	assert(p == N_fn);
    } 
    else if (case_number == 1) {
        int p = 0;

	label(p,0) = 1;    label(p++,1) = 0;
	label(p,0) = 1;    label(p++,1) = 1;
	label(p,0) = 1;    label(p++,1) = 2;

	label(p,0) = 1;    label(p++,1) = nbrs[1][N1-2];
	label(p,0) = 2;    label(p++,1) = nbrs[2][1];
	label(p,0) = 2;    label(p++,1) = 2;
	label(p,0) = 2;    label(p++,1) = 0;
	label(p,0) = 0;    label(p++,1) = 0;
	label(p,0) = 0;    label(p++,1) = nbrs[0][N0-2];
	for (int i = 1; i < N1 - 2; i++) {
	    label(p,0) = 1;    label(p++,1) = nbrs[1][i];
	}

	assert(p == N_fn);
    } 
    else if (case_number == 2) {
        int p = 0;

	label(p,0) = 2;    label(p++,1) = 0;
	label(p,0) = 2;    label(p++,1) = 1;
	label(p,0) = 2;    label(p++,1) = 2;

	for (int i = 1; i < N2 - 1; i++) {
	    label(p,0) = 2;    label(p++,1) = nbrs[2][i];
	}
	label(p,0) = 0;    label(p++,1) = nbrs[0][1];
	label(p,0) = 0;    label(p++,1) = 0;
	label(p,0) = 0;    label(p++,1) = 1;
	label(p,0) = 1;    label(p++,1) = 1;
	label(p,0) = 1;    label(p++,1) = nbrs[1][N1-2];

	assert(p == N_fn);
    } 
    else {
        int p = 0;

	label(p,0) = 1;    label(p++,1) = 2;
	label(p,0) = 2;    label(p++,1) = 0;
	label(p,0) = 0;    label(p++,1) = 1;

	label(p,0) = 0;    label(p++,1) = 0;
	label(p,0) = 0;    label(p++,1) = nbrs[0][N0-2];
	label(p,0) = 1;    label(p++,1) = nbrs[1][1];

	label(p,0) = 1;    label(p++,1) = 1;
	label(p,0) = 1;    label(p++,1) = nbrs[1][N1-2];
	label(p,0) = 2;    label(p++,1) = nbrs[2][1];

	label(p,0) = 2;    label(p++,1) = 2;
	label(p,0) = 2;    label(p++,1) = nbrs[2][N2-2];
	label(p,0) = 0;    label(p++,1) = nbrs[0][1];

	assert(p == N_fn);
    }

    // Distribute weights back to the coarse mesh
    // Init
    f.resize(N);
    df.resize(N,2);
    ddf.resize(N,2,2);

    f = 0.0;
    df = 0.0;
    ddf = 0.0;

    for (int p = 0; p < N_fn; p++) {
        int i0 = label(p,0);
	int i1 = label(p,1);

	if (i0 != i1) {
	    // The new point is an edge mid-point
	    int nnbr = nbrs[i0].size();

	    int j0, j1;
	    for (int k = 0; k < nnbr; k++) {
	        if ( nbrs[i0][k] == i1 ) {
		    int k0 = (k > 0 ? k - 1 : nnbr-1);
		    j0 = nbrs[i0][k0];

		    int k1 = (k+1)%nnbr;
		    j1 = nbrs[i0][k1];

		    break;
		}
	    }

	    f(i0) += 0.375*f_fn(p);
	    f(i1) += 0.375*f_fn(p);
	    f(j0) += 0.125*f_fn(p);
	    f(j1) += 0.125*f_fn(p);

	    m_daxpy(2, 0.375, &df_fn(p,0), &df(i0,0));
	    m_daxpy(2, 0.375, &df_fn(p,0), &df(i1,0));
	    m_daxpy(2, 0.125, &df_fn(p,0), &df(j0,0));
	    m_daxpy(2, 0.125, &df_fn(p,0), &df(j1,0));

	    m_daxpy(4, 0.375, &ddf_fn(p,0,0), &ddf(i0,0,0));
	    m_daxpy(4, 0.375, &ddf_fn(p,0,0), &ddf(i1,0,0));
	    m_daxpy(4, 0.125, &ddf_fn(p,0,0), &ddf(j0,0,0));
	    m_daxpy(4, 0.125, &ddf_fn(p,0,0), &ddf(j1,0,0));
	}
	else {
	    // The new point is still a vertex of the coarser mesh
	    int nnbr = nbrs[i0].size();
	    double wght = 1.0/nnbr*(0.625 - pow(0.375 + 0.25*cos(2*M_PI/nnbr), 2));
	    double wght_center = 1.0 - nnbr*wght;

	    for (int k = 0; k < nnbr; k++) {
		int i = nbrs[i0][k];
		f(i) += wght*f_fn(p);
		m_daxpy(2, wght, &df_fn(p,0), &df(i,0));
		m_daxpy(4, wght, &ddf_fn(p,0,0), &ddf(i,0,0));
	    }

	    f(i0) += wght_center*f_fn(p);
	    m_daxpy(2, wght_center, &df_fn(p,0), &df(i0,0));
	    m_daxpy(4, wght_center, &ddf_fn(p,0,0), &ddf(i0,0,0));
	}
    } // p

    // Scale the gradients
    if (case_number != 3) 
	df *= 2.0;
    else
	df *= -2.0;

    ddf *= 4.0;
}


/* Calculate shape function with the option of tabulating the results 
 * Arguments:
 *  N0,N1,N2 -- valance of the three vertices
 *  v,w -- barycentric coordinate
 *  f, df, ddf -- 
 *  updateDb -- whether to update the database with the result */
void subdiv::calcFunc(int N0, int N1, int N2, double v, double w, 
	MArray<double,1> &f, MArray<double,2> &df, MArray<double,3> &ddf,
	bool updateDb)
{
    static Db db;
    DbKey key(N0,N1,N2,v,w);
    Db::iterator it = db.find(key);

    int N = N0 + N1 + N2 - 6;
    f.resize(N);
    df.resize(N,2);
    ddf.resize(N,2,2);

    if (it != db.end()) {
        DbVal &val = (*it).second;

        std::copy(val.begin(), val.begin()+N, f.data());
	std::copy(val.begin()+N, val.begin()+3*N, df.data());
	std::copy(val.begin()+3*N, val.begin()+7*N, ddf.data());
    }
    else if (updateDb) {
        MArray<double,1> f_db(N);
        MArray<double,2> df_db(N,2);
        MArray<double,3> ddf_db(N,2,2);
	calcFunc_Recursive(N0, N1, N2, v, w, f_db, df_db, ddf_db);

	// Update database
	vector<double> val(7*N);
	std::copy(f_db.data(), f_db.data()+N, val.begin());
	std::copy(df_db.data(), df_db.data()+2*N, val.begin()+N);
	std::copy(ddf_db.data(), ddf_db.data()+4*N, val.begin()+3*N);
	db[key] = val;

	f = f_db;
	df = df_db;
	ddf = ddf_db;
    }
    else {
	calcFunc_Recursive(N0, N1, N2, v, w, f, df, ddf);
    }
}


/* Build subdivision function data base 
 * Arguments:
 *  qrule -- quadrature rule */
void subdiv::buildFuncDataBase(const Quad2D &qrule)
{
    MArray<double,1> f;
    MArray<double,2> df;
    MArray<double,3> ddf;

    for (int iq = 0; iq < qrule.n(); iq++) {

	for (int N0 = 4; N0 <= 8; N0++)
	for (int N1 = 4; N1 <= 8; N1++)
	for (int N2 = 4; N2 <= 8; N2++) {
	    calcFunc(N0, N1, N2, qrule.x(iq), qrule.y(iq), f, df, ddf, true);
	}
    }
}


/* For each triangle, find the lists of vertices that has non-zero shape 
 * function inside the triangle 
 * Arguments:
 *  nvert -- number of verts
 *  nface -- number of faces 
 *  ctrlPts -- all control points */
void subdiv::buildControlPoints(int nvert, int nface, const int (*f2v)[3], 
		vector<vector<int> > &ctrlPts)
{
    // Build the one-ring neighbor list of all vertices
    vector<vector<int> > ring;
    Mesh_vertOneRingNbrs(nvert, nface, f2v, ring);

    ctrlPts.resize(nface);

    for (int iface = 0; iface < nface; iface++) {
	vector<int> &list_tmp = ctrlPts[iface];
	list_tmp.clear();

	int i0 = f2v[iface][0];
	int i1 = f2v[iface][1];
	int i2 = f2v[iface][2];

	list_tmp.push_back(i0);
	list_tmp.push_back(i1);
	list_tmp.push_back(i2);

	vector<int> &ring0 = ring[i0];
	vector<int> &ring1 = ring[i1];
	vector<int> &ring2 = ring[i2];

	vector<int>::iterator p;
        p = std::find(ring2.begin(), ring2.end(), i1);
	assert( p != ring2.end() );
	std::rotate(ring2.begin(), p, ring2.end());
	p = std::find(ring2.begin(), ring2.end(), i0);
	assert( p != ring2.end() );
	list_tmp.insert(list_tmp.end(), ring2.begin()+1, p);

        p = std::find(ring0.begin(), ring0.end(), list_tmp.back());
	assert( p != ring0.end() );
	std::rotate(ring0.begin(), p, ring0.end());
	p = std::find(ring0.begin(), ring0.end(), i1);
	assert( p != ring0.end() );
	list_tmp.insert(list_tmp.end(), ring0.begin()+1, p);

        p = std::find(ring1.begin(), ring1.end(), list_tmp.back());
	assert( p != ring1.end() );
	std::rotate(ring1.begin(), p, ring1.end());
	// The first vertex on the outer ring is list_tmp[3]
	p = std::find(ring1.begin(), ring1.end(), list_tmp[3]);
	assert( p != ring1.end() );
	list_tmp.insert(list_tmp.end(), ring1.begin()+1, p);
    } // iface
}
