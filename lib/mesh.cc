#include "cxxheaders.h"
#include "mesh.h"
#include "ewald.h"
#include "mathfunc.h"
#include "hdf5.h"
#include "hdf5_hl.h"

/* Mesh::Mesh 
 * Mesh::Mesh(const Mesh&)
 * Mesh::operator=(const Mesh&)
 * Mesh::~Mesh 
 * Mesh::getCoords 
 * Mesh::setCoords 
 * Mesh::getConnectivities 
 * Mesh::setConnectivities 
 * Mesh::setInternalPointers 
 * Mesh::connectPeriodicBoundaries 
 * Mesh::noPeriodicBoundaries 
 * Mesh::updateGeometry
 * Mesh::updateAVC
 * Mesh::getCoordRange 
 * Mesh::calcVertArea 
 * Mesh::calcVertNormal
 * Mesh::calcMomentInertia
 * Mesh::shapeFactor
 * Mesh::calcCentroidVel
 * Mesh::calcTransRotatVel
 * Mesh::calcMinMaxEdgeLen
 * Mesh::minEdgeLen
 * Mesh::maxEdgeLen
 * Mesh::minAngle
 * Mesh::isInteriorPoint
 * Mesh::copyMeshFrom 
 * Mesh::buildEdgeList 
 * Mesh::checkOrientation
 * Mesh::calcDoubleLayerLinearTerm 
 * Mesh::calcDoubleLayerJump 
 * Mesh::tensionForce 
 * Mesh::velDiv 
 * Mesh::checkConnectivity 
 * Mesh::readHDF5
 * Mesh::readDATsingle
 * Mesh::writeHDF5
 * Mesh::writeTecplot
 * Mesh_edgeList 
 * Mesh_vertOneRingNbrs
 * Mesh_faceNbrList
 * Mesh_reorderF2v  
 * Mesh_slice */

extern "C" void DSYEV_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*);


/* Edge link list */
namespace {
    struct EdgeWithLnk {
	int _v0, _v1;	// vertex index
	int _next;	// next edge

	EdgeWithLnk() { _v0 = _v1 = _next = -1; }
	EdgeWithLnk(int v0, int v1) : _v0(v0), _v1(v1), _next(-1) { }
	EdgeWithLnk(int v0, int v1, int next) : _v0(v0), _v1(v1), _next(next) { }

	int v0() const { return _v0; }
	int v1() const { return _v1; }
	int next() const { return _next; }
    };
}

Mesh::Mesh()
{
}


Mesh::~Mesh()
{
}


Mesh::Mesh(const Mesh &mesh) 
{
    int nvert = mesh.numVerts();
    verts.resize(nvert);
    for (int ivert = 0; ivert < nvert; ivert++) {
	FOR_I3 verts[ivert].x[i] = mesh.verts[ivert].x[i];
    }

    int nface = mesh.numFaces();
    faces.resize(nface);
    for (int iface = 0; iface < nface; iface++) {
        for (int l = 0; l < 3; l++) 
	    faces[iface].ivert[l] = mesh.faces[iface].ivert[l];
    }
}


Mesh& Mesh::operator=(const Mesh &mesh)
{
    int nvert = mesh.numVerts();
    verts.resize(nvert);
    for (int ivert = 0; ivert < nvert; ivert++) {
	FOR_I3 verts[ivert].x[i] = mesh.verts[ivert].x[i];
    }

    int nface = mesh.numFaces();
    faces.resize(nface);
    for (int iface = 0; iface < nface; iface++) {
        for (int l = 0; l < 3; l++) 
	    faces[iface].ivert[l] = mesh.faces[iface].ivert[l];
    }

    return *this;
}


/* Get all point coordinates
 * Arguments:
 *   x -- vertex coordinates stored in 1-d format */
void Mesh::getCoords(double *x) const
{
    int nvert = numVerts();

    int p = 0;
    for (int ivert = 0; ivert < nvert; ivert++) {
        const Point &vert = verts[ivert];
	FOR_I3 x[p++] = vert.x[i];
    }
}

void Mesh::getCoords(double (*x)[3]) const
{
    getCoords(*x);
}

void Mesh::getCoords(MArray<double,2> &x) const
{
    x.resize(numVerts(), 3);
    getCoords(x.data());
}


/* Set mesh point coordinates
 * Arguments:
 *   x -- vertex coordinate in 1D format*/
void Mesh::setCoords(const double *x)
{
    int nvert = numVerts();

    int p = 0;
    for (int ivert = 0; ivert < nvert; ivert++) {
        Point &vert = verts[ivert];
	FOR_I3 vert.x[i] = x[p++];
    }
}


void Mesh::setCoords(const double (*x)[3])
{
    setCoords(*x);
}

void Mesh::setCoords(const MArray<double,2> &x)
{
    setCoords(x.data());
}


/* Get mesh connectivities
 * Arguments:
 *   f2v -- face to vertex connectivities */
void Mesh::getConnectivities(int *f2v) const
{
    int nface = numFaces();

    int p = 0;
    for (int iface = 0; iface < nface; iface++) {
        const Tri &face = faces[iface];
	FOR_I3 f2v[p++] = face.ivert[i];
    }
}


void Mesh::getConnectivities(int (*f2v)[3]) const
{
    getConnectivities(*f2v);
}


void Mesh::getConnectivities(MArray<int,2> &f2v) const
{
    f2v.resize(numFaces(),3);
    getConnectivities(f2v.data());
}


/* Set mesh connectivity
 * Arguments:
 *   f2v -- face to vertex connectivity table */
void Mesh::setConnectivities(const int *f2v)
{
    int nface = numFaces();

    int p = 0;
    for (int iface = 0; iface < nface; iface++) {
        Tri &face = faces[iface];
	FOR_I3 face.ivert[i] = f2v[p++];
    }
}


void Mesh::setConnectivities(const int (*f2v)[3])
{
    setConnectivities(*f2v);
}


void Mesh::setConnectivities(const MArray<int,2> &f2v)
{
    setConnectivities(f2v.data());
}


/* Set internal pointers
 *  -- vertex-to-mesh, face-to-mesh, face-to-vert */
void Mesh::setInternalPointers()
{
    // Vertex
    for (int ivert = 0; ivert < numVerts(); ++ivert) {
        Point &vert = verts[ivert];

        vert.mesh = this;
	vert.indx = ivert;
    }

    // Face
    for (int iface = 0; iface < numFaces(); ++iface) {
        Tri &face = faces[iface];

        face.mesh = this;
	FOR_I3 face.vert[i] = &verts[face.ivert[i]];
    }

    // Edge
    buildEdgeList();
}


/* Connect the points on periodic boundaries
 * Count the number of free mesh points
 * Arguments:
 *   L -- periodic domain size */
void Mesh::connectPeriodicBoundaries(const double *L)
{
    // Init
    for (int ivert = 0; ivert < numVerts(); ivert++) {
        verts[ivert].pp = NULL;
    }

    if (L == NULL) return;

    // Find mesh boundary
    double xmin[3], xmax[3], eps[3];
    bool is_periodic_bd[3];

    getCoordRange(xmin, xmax);
    FOR_I3 {
	eps[i] = 1.E-10*max(L[i], xmax[i]-xmin[i]);
	is_periodic_bd[i] = fabs(xmax[i] - xmin[i] - L[i]) < eps[i];
    }

    // Find all the mesh points on the boundary
    vector<Point*> vbd;

    for (int ivert = 0; ivert < numVerts(); ivert++) {
        Point &vert = verts[ivert];
	bool onBd = false;

	FOR_I3 {
	    if (! is_periodic_bd[i]) continue;
	    onBd = onBd || (vert.x[i] < xmin[i] + eps[i]);
	    onBd = onBd || (vert.x[i] > xmax[i] - eps[i]);
	}

	if (onBd) vbd.push_back(&vert);
    }

    // Connect points on the periodic boundaries
    for (int i0 = 0; i0 < vbd.size()-1; i0++) {
        Point *v0 = vbd[i0];

	for (int i1 = i0+1; i1 < vbd.size(); i1++) {
	    Point *v1 = vbd[i1];

	    bool same_point = true;

	    FOR_J3 {
		double xx = v1->x[j] - v0->x[j];
		xx = xx - nearbyint(xx/L[j])*L[j];
		if (fabs(xx) > eps[j]) {
		    same_point = false;
		    break;
		}
	    }

	    if (same_point) {
	        v1->pp = v0;
	        break;
	    }
	} // i1
    } // i0
}


/* No periodic boundaries */
void Mesh::noPeriodicBoundaries()
{
    for (int ivert = 0; ivert < numVerts(); ivert++) {
        verts[ivert].pp = NULL;
    }
}


/* Fully update geometry information */
void Mesh::updateGeometry()
{
    // Faces
    for (int iface = 0; iface < numFaces(); iface++) {
        faces[iface].updateGeometry();
    }
    
    // Area, centroid and volume
    area = 0.0;
    vol = 0.0;
    FOR_D3 center[d] = 0.0;

    // A reference point
    const double *x0 = verts[0].x;

    // Area and center
    for (int iface = 0; iface < numFaces(); ++iface) {
        Tri &face = faces[iface];
        area += face.area;

	double xx[3];
	FOR_I3 xx[i] = face.xc[i] - x0[i];
	vol += THRD*m_ddot(3, xx, face.normal)*face.area;

	FOR_I3 center[i] += face.area*face.xc[i];
    }
    FOR_I3 center[i] /= area;

    // Moment of inertia
    calcMomentInertia(pmom, paxis);

    // Vertex area, normal and curvature
    calcVertArea(vertArea);
    calcVertNormal(vertNrml, vertH);
}


/* Update the area, vol, and center 
 * Note:
 *   -- This function does not compute the detailed geometries of each
 *      face, so is more economical. */
void Mesh::updateAVC()
{
    // Init
    area = 0.0;
    vol = 0.0;
    FOR_D3 center[d] = 0.0;;

    // Reference point
    const double *x0 = verts[0].x;

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &T = faces[iface];

        double xtri[3][3]; 
	double xc[3] = {0.0};
	for (int l = 0; l < 3; l++) {
	    int ivert = T.ivert[l];
	    Point &vert = verts[ivert];

	    FOR_D3 xtri[l][d] = vert.x[d] - x0[d];
	    FOR_D3 xc[d] += THRD*vert.x[d];
	}

	double normal[3], detj, dA;
	tri_normal(xtri, normal, detj);
	dA = 0.5*detj;

	area += dA;
	vol += THRD*dA*m_ddot(3, normal, xtri[0]);
	FOR_D3 center[d] += dA*xc[d];
    }

    FOR_D3 center[d] /= area;
}


/* Calculate the coordinate range
 * Arguments:
 *   xmin, xmax -- minimum and maximum of the mesh coordinates */
void Mesh::getCoordRange(double *xmin, double *xmax)
{
    xmin[0] = xmin[1] = xmin[2] = FLT_MAX;
    xmax[0] = xmax[1] = xmax[2] = -FLT_MAX;

    for (int i = 0; i < numVerts(); ++i) {
        const double *xtmp = verts[i].x;

        for (int d = 0; d < 3; ++d) {
	    xmin[d] = min(xmin[d], xtmp[d]);
	    xmax[d] = max(xmax[d], xtmp[d]);
        } // d
    } // i
}


/* Calculate the coordinate range of a specific dimension
 * Arguments:
 *   d -- the dimension
 *   xmin, xmax */
void Mesh::getCoordRange(int d, double &xmin, double &xmax)
{
    assert(d >= 0 && d < 3);

    xmin = FLT_MAX;
    xmax = -FLT_MAX;

    for (int i = 0; i < numVerts(); ++i) {
        const double xtmp = verts[i].x[d];
	xmin = min(xmin, xtmp);
	xmax = max(xmax, xtmp);
    } // i
}


/* Compute the area of the surface occupied by each vertex
 * Arguments:
 *   va -- vertex area
 * Note:
 *   The vertex area is defined as 1/3 of the sum of the area of
 *   the elements surrounding a mesh vertex */
void Mesh::calcVertArea(MArray<double,1> &va)
{
    int nvert = numVerts();

    // Init
    va.resize(nvert);
    va = 0.0;

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &face = faces[iface];
	for (int l = 0; l < 3; l++) {
	    int ivert = face.ivert[l];
	    va(ivert) += face.area;
	}
    }

    va *= THRD;
}


/* Calculate the surface normal and mean curvature at each vertex position
 * Arguments:
 *   nrml -- surface normal
 *   H -- mean curvature */
void Mesh::calcVertNormal(MArray<double,2> &nrml, MArray<double,1> &H)
{
    int nvert = numVerts();
    int nface = numFaces();

    // Calc vertex area
    MArray<double,1> va(nvert);
    MArray<double,2> nrml_approx(nvert,3);		
    va = 0.0;
    nrml_approx = 0.0;

    for (int iface = 0; iface < nface; iface++) {
        Tri &face = faces[iface];
	for (int l = 0; l < 3; l++) {
	    int ivert = face.ivert[l];
	    va(ivert) += THRD*face.area;

	    m_daxpy(3, face.area, face.normal, &nrml_approx(ivert,0));
	}
    }

    for (int ivert = 0; ivert < nvert; ivert++) {
        normalizeVec3D( &nrml_approx(ivert,0) );
    }

    // Init
    nrml.resize(nvert,3);
    H.resize(nface);

    nrml = 0.0;
    H = 0.0;

    for (int iface = 0; iface < nface; iface++) {
        Tri &face = faces[iface];

	for (int l = 0; l < 3; l++) {
	    int iv0 = face.ivert[l];
	    int iv1 = face.ivert[(l+1)%3];
	    int iv2 = face.ivert[(l+2)%3];

	    double t[3];
	    FOR_D3 t[d] = 0.5*(verts[iv2].x[d] - verts[iv1].x[d]);

	    double txn[3];
	    cross_product(t, face.normal, txn);

	    FOR_D3 nrml(iv0,d) += txn[d];
	}
    }

    for (int ivert = 0; ivert < nvert; ivert++) {
        H(ivert) = normalizeVec3D( &nrml(ivert,0)  );
	H(ivert) /= 2*va(ivert);

	if ( m_ddot(3, &nrml(ivert,0), &nrml_approx(ivert,0)) < 0) {
	    FOR_D3 nrml(ivert,d) = -nrml(ivert,d);
	    H(ivert) = -H(ivert);
	}
    }
}


/* Compute the moment of inertia tensor 
 * Argument:
 *  lbd -- eigenvalues
 *  axis[i][:] -- i-th eigenvector
 * Note:
 *   Use the center of the mesh as origin */
void Mesh::calcMomentInertia(double *lbd, double (*axis)[3])
{
    // Calculate the moment of inertia tensor
    double mm[3][3] = { 0.0 };

    Quad2D &Q = quadrature::select_rule_2d("TRI_3");
    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &face = faces[iface];

	for (int iq = 0; iq < Q.n(); ++iq) {
	    double s, t, w0, w1, w2, dA;
	    double rr2, xx[3];

	    s = Q.x(iq);
	    t = Q.y(iq);

	    w0 = 1.0 - s - t;
	    w1 = s;
	    w2 = t;

	    dA = Q.w(iq)*face.detJ;

	    FOR_I3 xx[i] = w0*face.vert[0]->x[i] 
	                 + w1*face.vert[1]->x[i] 
			 + w2*face.vert[2]->x[i] 
			 - center[i];

	    rr2 = xx[0]*xx[0] + xx[1]*xx[1] + xx[2]*xx[2];

	    FOR_I3 {
	        mm[i][i] += dA*rr2;
	        FOR_J3 mm[i][j] -= dA*xx[i]*xx[j];
	    }
	} // iq
    } // iface

    // Force symmetry
    mm[1][0] = mm[0][1];
    mm[2][0] = mm[0][2];
    mm[2][1] = mm[1][2];

    // Find the principal axes and eigenvalues
    char JOBZ = 'V';
    char UPLO = 'U';
    int N = 3;
    double A[9];
    int LDA = 3;
    double W[3];
    int LWORK = 15;
    double WORK[15];
    int INFO;

    // Eigen-analysis by Lapack subroutine DSYEV_
    // Copy mm(:,:) to A (Fortran style array)
    FOR_I3 FOR_J3 A[i+3*j] = mm[i][j];
    dsyev_(&JOBZ, &UPLO, &N, A, &LDA, W, WORK, &LWORK, &INFO);

    FOR_I3 lbd[i] = W[i];
    // Eigenvectors are saved as columns in A (Fortran style array)
    // Hence one vector is stored after another in memory continuously 
    //
    // axis[i][:] = A(:,i) (right-hand side Fortran array)
    m_dcopy(9, A, *axis); 	
}


/* Compute the shape factor
 * Arguments:
 *  psi -- inclination angle
 *  D -- Taylor deformation factor
 * Note:
 *  Assume that the x-z plane is the shear plane */
void Mesh::shapeFactor(double &psi, double &D)
{
    // Length of the vesicle along the principle axes
    double Lmin[3], Lmax[3], L[3];
    FOR_D3 {
        Lmin[d] = FLT_MAX;
        Lmax[d] = -FLT_MAX;
    }

    for (int ivert = 0; ivert < numVerts(); ivert++) {
        Point &vert = verts[ivert];
	FOR_I3 {
	    double Ltmp = m_ddot(3, vert.x, paxis[i]);
	    Lmin[i] = min(Lmin[i], Ltmp);
	    Lmax[i] = max(Lmax[i], Ltmp);
	}
    }

    FOR_I3 L[i] = Lmax[i] - Lmin[i];

    // Indentify the major axis
    int ix = 0;
    for(int i = 1; i < 3; i++) {
	if (L[i] > L[ix]) ix = i;
    }

    // Identify the axis that is parallel to the y-axis (i.e. the vorticity direction)
    // The third axis is the z-axis
    int iy = (ix + 1)%3;
    int iz = (iy + 1)%3;
    if ( fabs(paxis[iz][1]) > fabs(paxis[iy][1]) ) swap(iy, iz);


    psi = atan2(paxis[ix][2], paxis[ix][0]);
    psi -= nearbyint(psi/M_PI)*M_PI;

    D = (L[ix] - L[iz])/(L[ix] + L[iz]);
}


/* Compute the velocity of the volumetric center
 * Arguments:
 *   u -- the surface velocity
 *   UT -- the centroid velocity 
 * Note:
 *   -- Assume that the velocity field is divergence free
 *  */
void Mesh::calcCentroidVel(const double *u, double *UT)
{
    // Init
    FOR_D3 UT[d] = 0.0;

    Quad2D &Q = quadrature::select_rule_2d("TRI_3");

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &tri = faces[iface];
        double xtri[3][3], utri[3][3];
	double xq[3], uq[3];

        for (int l = 0; l < 3; l++) {
            int &ivert = tri.ivert[l];
	    m_dcopy(3, verts[ivert].x, xtri[l]);
	    m_dcopy(3, u+3*ivert, utri[l]);
        }

        for (int iq = 0; iq < Q.n(); iq++) {
            double s=Q.x(iq), t=Q.y(iq);
	    double r=1.0-s-t;
	    double dA=tri.detJ*Q.w(iq);
            double xx[3];

            FOR_I3 {
                xx[i] = r*xtri[0][i] + s*xtri[1][i] + t*xtri[2][i] - center[i];
                uq[i] = r*utri[0][i] + s*utri[1][i] + t*utri[2][i];
            }

	    double uqn = m_ddot(3, uq, tri.normal);
	    m_daxpy(3, uqn*dA, xx, UT);
        }
    } // iface

    // Normalize
    FOR_I3 UT[i] /= vol;
}


void Mesh::calcCentroidVel(const double (*u)[3], double *UT)
{
    calcCentroidVel(*u, UT);
}


void Mesh::calcCentroidVel(const MArray<double,2> &u, double *UT)
{
    calcCentroidVel(u.data(), UT);
}


/* Compute the translational and rotational velocity
 * Arguments:
 *   u -- the surface velocity
 *   UT -- translating velocity
 *   OMG -- angular velocity 
 * Note:
 *   -- The principle axis paxis[][] and principle momentum pmm[] must 
 *      first be updated */
void Mesh::calcTransRotatVel(const double *u, double *UT, double *OMG)
{
    Quad2D &Q = quadrature::select_rule_2d("TRI_3");

    double UR[3];	// angular velocity along each principle axis

    // Init
    FOR_D3 {
        UT[d] = 0.0;
	UR[d] = 0.0;
    }

    // Calculate inner-product with rigid body modes
    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &tri = faces[iface];
        double xtri[3][3], utri[3][3];
	double xq[3], uq[3];

        for (int l = 0; l < 3; l++) {
            int &ivert = tri.ivert[l];

	    m_dcopy(3, verts[ivert].x, xtri[l]);
	    m_dcopy(3, u+3*ivert, utri[l]);
        }

        for (int iq = 0; iq < Q.n(); iq++) {
            double s=Q.x(iq), t=Q.y(iq), dA=tri.detJ*Q.w(iq);
            double w0=1.0-s-t, w1=s, w2=t;
            double xx[3], gg[3], ggnn;

            FOR_I3 {
                xx[i] = w0*xtri[0][i] + w1*xtri[1][i] + w2*xtri[2][i] - center[i];
                uq[i] = w0*utri[0][i] + w1*utri[1][i] + w2*utri[2][i];
            }

	    // Translation
	    FOR_I3 UT[i] += uq[i]*dA;

	    // Rotation
	    FOR_I3 {
	        double omgxx[3];
		cross_product(paxis[i], xx, omgxx);
		UR[i] += m_ddot(3, omgxx, uq)*dA;
	    }
        }
    } // iface

    // Normalize
    FOR_I3 UT[i] /= area;
    FOR_I3 UR[i] /= pmom[i];

    // Convert angular velocity to Cartesian coordinates
    FOR_I3 OMG[i] = UR[0]*paxis[0][i] + UR[1]*paxis[1][i] + UR[2]*paxis[2][i];
}


void Mesh::calcTransRotatVel(const double (*u)[3], double *UT, double *OMG)
{
    calcTransRotatVel(*u, UT, OMG);
}


void Mesh::calcTransRotatVel(const MArray<double,2> &u, double *UT, double *OMG)
{
    calcTransRotatVel(u.data(), UT, OMG);
}


void Mesh::calcMinMaxEdgeLen(double *minlen, double *maxlen)
{
    if (minlen) *minlen = FLT_MAX;
    if (maxlen) *maxlen = 0.0;

    int nface = numFaces();
    for (int iface = 0; iface < nface; iface++) {
        Tri &face = faces[iface];

        for (int l = 0; l < 3; l++) {
	    int iv0 = face.ivert[l];
	    int iv1 = face.ivert[(l+1)%3];
	    if (iv0 > iv1) continue;

	    double len = distCart3D(verts[iv0].x, verts[iv1].x);

	    if (minlen) *minlen = min(*minlen, len);
	    if (maxlen) *maxlen = max(*maxlen, len);
	}
    } //iface
}


/* The smallest edge length */
double Mesh::minEdgeLen()
{
    double minlen;
    calcMinMaxEdgeLen(&minlen, NULL);
    return minlen;
}


/* The max edge length */
double Mesh::maxEdgeLen()
{
    double maxlen;
    calcMinMaxEdgeLen(NULL, &maxlen);
    return maxlen;
}


/* The smallest angle */
double Mesh::minAngle()
{
    double thmin = FLT_MAX;	// initialize

    for (int ifa = 0; ifa < numFaces(); ifa++) {
        Tri &face = faces[ifa];

	double x01[3], x02[3], x12[3];
	for (int ii = 0; ii < 3; ii++) {
	    x01[ii] = face.vert[1]->x[ii] - face.vert[0]->x[ii];
	    x02[ii] = face.vert[2]->x[ii] - face.vert[0]->x[ii];
	    x12[ii] = face.vert[2]->x[ii] - face.vert[1]->x[ii];
	}

	double l01, l02, l12;
	l01 = m_dnrm2(3, x01);
	l02 = m_dnrm2(3, x02);
	l12 = m_dnrm2(3, x12);

	double th0, th1, th2;
	th0 = m_ddot(3, x01, x02)/(l01*l02);
	th0 = acos(th0);

	th1 = -m_ddot(3, x01, x12)/(l01*l12);
	th1 = acos(th1);

	th2 = M_PI - th0 - th1;

	thmin = min(thmin, th0);
	thmin = min(thmin, th1);
	thmin = min(thmin, th2);
    }

    return thmin;
}


/* Tell if a point is the interior of the mesh
 * Arguments:
 *   x -- the point */
bool Mesh::isInteriorPoint(const double *x)
{
    double xmin[3], xmax[3];
    getCoordRange(xmin, xmax);
    if (   x[0] < xmin[0] || x[0] > xmax[0] 
	|| x[1] < xmin[1] || x[1] > xmax[1]
	|| x[2] < xmin[2] || x[2] > xmax[2] ) return false;

    double sangle = 0.0;
    for (int iface = 0; iface < numFaces(); iface++) {
	Tri &face = faces[iface];
	double xtri[3][3];
	FOR_I3 FOR_J3 xtri[i][j] = face.vert[i]->x[j];
	sangle += tri_solidAngle(x, xtri);
    }

    if (fabs(sangle) > 1.E-2) 
	return true;
    else
	return false;
}


// Copy the mesh to another mesh
// Arguments:
//   smesh -- the source mesh
void Mesh::copyMeshFrom(const Mesh &smesh)
{
    int nvert = smesh.numVerts();
    int nface = smesh.numFaces();

    verts.resize(nvert);
    faces.resize(nface);

    for (int ivert = 0; ivert < nvert; ivert++) {
        FOR_I3 verts[ivert].x[i] = smesh.verts[ivert].x[i];
    }

    for (int iface = 0; iface < nface; iface++) {
        FOR_I3 faces[iface].ivert[i] = smesh.faces[iface].ivert[i];
    }
}


// Build edge list
void Mesh::buildEdgeList()
{
    int nvert = numVerts();
    int nface = numFaces();
    int nedgeMax = 3*nface;
    int nedge;

    // Allocate working arrays
    int (*f2v)[3] = new int[nface][3];
    int (*e2v)[2] = new int[nedgeMax][2];
    int (*e2f)[2] = new int[nedgeMax][2];

    // Build edge list
    getConnectivities(f2v);
    Mesh_edgeList(nvert, nface, f2v, nedge, e2v, e2f);

    // Copy edge list
    edges.resize(nedge);
    for (int iedge = 0; iedge < nedge; ++iedge) {
        Edge &edge = edges[iedge];

	for (int l = 0; l < 2; ++l) {
	    int ivert = e2v[iedge][l];
            edge.ivert[l] = ivert;
	    edge.vert[l] = (ivert >= 0)? &verts[ivert] : NULL;

	    int iface = e2f[iedge][l];
            edge.iface[l] = iface;
	    edge.face[l] = (iface >= 0)? &faces[iface] : NULL;
	}
    }

    // Deallocate working arrays
    delete []f2v;
    delete []e2v;
    delete []e2f;
}


void Mesh::checkOrientation()
{
    int nvert = numVerts();
    int nface = numFaces();
    int nedge;
    int NEDGE_MAX = 3*nface;
    int (*f2v)[3] = new int[nface][3];
    int (*e2v)[2] = new int[NEDGE_MAX][2];
    int (*e2f)[2] = new int[NEDGE_MAX][2];

    printf("Check mesh orientation\n");

    getConnectivities(f2v);
    Mesh_edgeList(nvert, nface, f2v, nedge, e2v, e2f);

    for (int iedge = 0; iedge < nedge; iedge++) {
        if (e2f[iedge][0] < 0 || e2f[iedge][1] < 0) continue;

        // We have make each edge to have to same orientation as its
	// first face neighbor.  It hence must have opposite orientation
	// as its second face neighbor
        bool pass0 = false;
	bool pass1 = false;

	int iface = e2f[iedge][0];
	for (int l = 0; l < 3; l++) {
	    if (f2v[iface][l] == e2v[iedge][0] && f2v[iface][(l+1)%3] == e2v[iedge][1]) {
		pass0 = true;
		break;
	    }
	}

	iface = e2f[iedge][1];
	for (int l = 0; l < 3; l++) {
	    if (f2v[iface][(l+1)%3] == e2v[iedge][0] && f2v[iface][l] == e2v[iedge][1]) {
		pass1 = true;
		break;
	    }
	}

	if ( !pass0 || !pass1 ) {
	    int iface0 = e2f[iedge][0];
	    int iface1 = e2f[iedge][1];
	    printf("  Faces %d(%d,%d,%d) and %d(%d,%d,%d) incompatible\n",
	    		iface0, f2v[iface0][0], f2v[iface0][1], f2v[iface0][2],
	    		iface1, f2v[iface1][0], f2v[iface1][1], f2v[iface1][2]);
	}
    }

    delete [] f2v;
    delete [] e2v;
    delete [] e2f;
}


/*
void Mesh::checkConnectivity()
{
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    if (edges.size() == 0) {
        if (mpi_rank == 0) {
	    printf("Edge list is empty.\n");
	    printf("No check performed.\n");
        }
	return;
    }

    // Build face to edge list
    int nface = numFaces();
    int (*f2e)[3] = new int[nface][3];
    for (int iface = 0; iface < nface; iface++) {
        f2e[iface][0] = -1;
        f2e[iface][1] = -1;
        f2e[iface][2] = -1;
    }

    for (int ied = 0; ied < edges.size(); ied++) {
	Edge &edge = edges[ied];
        int iv0 = edge.ivert[0];
	int iv1 = edge.ivert[1];

	for (int l = 0; l < 2; l++) {
	    int iface = edge.iface[l];
	    Tri &face = faces[iface];

	    bool flag = false;
	    for (int l0 = 0; l0 < 3; l0++) {
		int l1 = (l0+1)%3;
		if ( (face.ivert[l0] == iv0 && face.ivert[l1] == iv1) ||
		     (face.ivert[l0] == iv1 && face.ivert[l1] == iv0) ) {
		    f2e[iface][l0] = ied;
		    flag = true;
		    break;
		}
	    }
	    assert(flag);
        }
    }

    for (int iface = 0; iface < nface; iface++) {
        if (f2e[iface][0] < 0 || f2e[iface][1] < 0 || f2e[iface][2] < 0 ) {
	    printf("%d th proc\n", mpi_rank);
	    printf("Error: %d th face doesn't have three neighbors\n", iface);
	    throw(-1);
	}
    } //iface

    delete [] f2e;
}
*/


/* Calculate the linear term in double-layer operator
 * Arguments:
 *  rhs -- the linear term
 * Note:
 *  1. The double layer density is the g[] array of each vetex 
 */
void Mesh::calcDoubleLayerLinearTerm(double *rhs)
{
    Quad2D &Q = quadrature::select_rule_2d("TRI_3");
    
    // Init
    rhs[0] = rhs[1] = rhs[2] = 0.0;

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &tri = faces[iface];
	double xtri[3][3], gtri[3][3];

	for (int l = 0; l < 3; l++) {
	    int ivert = tri.ivert[l];
	    Point &vert = verts[ivert];

	    m_dcopy(3, vert.x, xtri[l]);
	    m_dcopy(3, vert.g, gtri[l]);
	}

	for (int iq = 0; iq < Q.n(); iq++) {
	    double s = Q.x(iq), t = Q.y(iq), dA = tri.detJ*Q.w(iq); 
	    double w0 = 1.0-s-t, w1 = s, w2 = t;
	    double xx[3], gq[3], gn;

	    FOR_I3 {
	        xx[i] = w0*xtri[0][i] + w1*xtri[1][i] + w2*xtri[2][i] - center[i];
		gq[i] = w0*gtri[0][i] + w1*gtri[1][i] + w2*gtri[2][i];
	    }

	    gn = m_ddot(3, gq, tri.normal);
	    m_daxpy(3, dA*gn, xx, rhs);
	}
    }

    m_dscal(3, -8.0*M_PI*ewald::iVol, rhs);
}


/* Calculate the double-layer princple value tensor
 * Note:
 *   1. For x0 on D, lim x->x0+ K*f(x) = 4*PI*C * f(x0) + K*f(x0)
 *      where K is the double layer operator
 *   2. For a smooth operator, C equals to the identity operator
 *   3. C is stored at every vertex */
void Mesh::calcDoubleLayerJump()
{
    assert(edges.size() > 0);

    int nvert = numVerts();
    double *SA = new double[nvert];	// SA = inner solid angles

    // Init
    for (int ivert = 0; ivert < nvert; ivert++) SA[ivert] = 2.0*M_PI;

    for (int iedge = 0; iedge < numEdges(); iedge++) {
	Edge &edge = edges[iedge];
	int ivert0 = edge.ivert[0];
	int ivert1 = edge.ivert[1];
	int iface0 = edge.iface[0];
	int iface1 = edge.iface[1];

        double axis[3], rr;
	FOR_I3 axis[i] = verts[ivert1].x[i] - verts[ivert0].x[i];
	normalizeVec3D(axis);

	double *normal0 = faces[iface0].normal;
	double *normal1 = faces[iface1].normal;
	double costh = m_ddot(3, normal0, normal1);
	double n01[3];
	cross_product(normal0, normal1, n01);
	double sinth = m_ddot(3, n01, axis);
	double dth = atan2(sinth, costh);

	SA[ivert0] -= dth;
	SA[ivert1] -= dth;
    }


    //======================================================================
    // Now compute the Cij tensor
    // Init
    for (int ivert = 0; ivert < numVerts(); ivert++) {
        Point &vert = verts[ivert];
        m_dclear(9, *vert.C);
	FOR_I3 vert.C[i][i] = SA[ivert]/(2*M_PI);
    }

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &face = faces[iface];
	for (int l = 0; l < 3; l++) {
	    int ivert0 = face.ivert[l];
	    int ivert1 = face.ivert[(l+1)%3];
	    int ivert2 = face.ivert[(l+2)%3];

	    Point &vert0 = verts[ivert0];
	    Point &vert1 = verts[ivert1];
	    Point &vert2 = verts[ivert2];

	    double axis1[3], axis2[3], rr;
	    FOR_I3 {
	        axis1[i] = vert1.x[i] - vert0.x[i];
	        axis2[i] = vert2.x[i] - vert0.x[i];
	    }

	    normalizeVec3D(axis1);
	    normalizeVec3D(axis2);

	    double th;
	    th = m_ddot(3, axis1, axis2);
	    th = min(1.0, max(-1.0, th));
	    th = acos(th);

	    double xc[3];
	    FOR_I3 xc[i] = 0.5*(axis1[i] + axis2[i]);
	    rr = m_dnrm2(3, xc);
            m_dscal(3, (1.0/M_PI)*sin(0.5*th)/rr, xc);

	    // Add to the Cij tensor at ivert0-th vertex
	    FOR_I3
	    FOR_J3 {
	        vert0.C[i][j] -= face.normal[i]*xc[j];
	    }
	} // l
    } // face

    // Dealloc temp arrays
    delete [] SA;
}


/* Compute the surface tension force
 * Arguments:
 *  sigma -- surface tension
 *  f -- tension force */
void Mesh::tensionForce(const MArray<double,1> &sigma, MArray<double,2> &f)
{
    int nvert = numVerts();

    // Init
    f.resize(nvert,3);
    f = 0.0;

    // f_i = 1/A_i * dA/dx
    for (int iface = 0; iface < numFaces(); ++iface) {
        Tri &T = faces[iface];

        double dAdx[3][3];
	T.gradArea(dAdx);

	double coeff = 0.0;
	for (int l = 0; l < 3; l++) {
	    int ivert = T.ivert[l];
	    coeff += THRD*sigma(ivert);
	}

	for (int l = 0; l < 3; l++) {
	    int ivert = T.ivert[l];
	    m_daxpy(3, coeff, dAdx[l], &f(ivert,0));
	}
    } //iface

    for (int ivert = 0; ivert < nvert; ivert++) {
        m_dscal(3, 1.0/vertArea(ivert), &f(ivert,0));
    }
}


/* Compute the surface divergence of the surface velocity field
 * Arguments:
 *  v -- the surface velocity 
 *  div -- velocity divergence */
void Mesh::velDiv(const MArray<double,2> &v, MArray<double,1> &div)
{
    // Init
    div = 0.0;

    for (int iface = 0; iface < numFaces(); iface++) {
        Tri &T = faces[iface];

        double dAdx[3][3];
        T.gradArea(dAdx);

        double tmp = 0.0;
        for (int l = 0; l < 3; l++) {
            int ivert = T.ivert[l];
            tmp += m_ddot(3, dAdx[l], &v(ivert,0));
        }
        tmp *= THRD;

        for (int l = 0; l < 3; l++) {
            int ivert = T.ivert[l];
            div(ivert) += tmp;
        }
    } // iface
}




/* Read HDF5 mesh */
void Mesh::readHDF5(const char *fn)
{
    hsize_t dims[2]; 
    H5T_class_t class_id;
    size_t type_size;
    int nvert, nface;
    double *xbuf;
    int *f2vbuf;

    hid_t fid = H5Fopen(fn, H5F_ACC_RDONLY, H5P_DEFAULT);

    H5LTget_dataset_info(fid, "X", dims, &class_id, &type_size);
    nvert = dims[0];

    H5LTget_dataset_info(fid, "F2V", dims, &class_id, &type_size);
    nface = dims[0];

    // Alloc temp arrays
    xbuf = new double[3*nvert];
    f2vbuf = new int[3*nface];

    H5LTread_dataset_double(fid, "X", xbuf);
    H5LTread_dataset_int(fid, "F2V", f2vbuf);

    verts.resize(nvert);
    for (int ivert = 0; ivert < nvert; ivert++) {
        Point &vert = verts[ivert];
	FOR_J3 vert.x[j] = xbuf[3*ivert+j];
    }

    faces.resize(nface);
    for (int iface = 0; iface < nface; iface++) {
        Tri &face = faces[iface];
	FOR_J3 face.ivert[j] = f2vbuf[3*iface+j];
    }

    // Delete temp arrays
    delete [] xbuf;
    delete [] f2vbuf;

    H5Fclose(fid);
}


//Function added by spann October 2012, analogous to readDAT in spann's singleVesicle code branch
/* Read .dat tecplot compatible mesh, such as those generated by singleVesicle code.  The file does not need to have the velocity information to be a valid input file; this routine will ignore any velocities that are present. */
void Mesh::readDATsingle(const char *fn)
{
    int nvert, nface;
    double *x;
    int *f2v;
	int pos = 0;
	string sline;
	ifstream infile;
	istringstream mystringstream;


	infile.open(fn);
	getline(infile,sline);   //Throw away 1st line, which is "# TIME ="
	getline(infile,sline);   //Throw away 2nd line, which reads XYZ (optional UVW)
	getline(infile,sline);   //Read third line to find nvert and nface.  begins with "ZONE N="
	sline = sline.substr(7); //Discard "ZONE N="
	mystringstream.clear();
	mystringstream.str(sline);  //Read the number of vertices
	mystringstream >> nvert;  //Read nvert
	//cout << "nvert = " << nvert << endl;
	sline = sline.substr(int(floor(log10(nvert)))+4);//Jump to number past E=
	mystringstream.clear();
	mystringstream.str(sline);  //Read the number of faces
	mystringstream >> nface;
    //cout << "nface = " << nface << endl;

	//nface = 2*nvert-4;        //nface computed from Euler's formula if all faces are triangles.  //Update:  not assuming this anymore because we could want to load walls


    x = new double[3*nvert];  //We're going to fill these arrays then transfer the values to the mesh's arrays
    f2v = new int[3*nface];

    	pos = 0;
		for(int i = 0; i < nvert; ++i)  //Read x coordinates, ignore velocities if present in the file
		{
			getline(infile,sline);
			mystringstream.clear();
			mystringstream.str(sline);
			mystringstream >> x[pos] >> x[pos+1] >> x[pos+2];  //first three values of each line are XYZ position
			pos += 3;
		}

		pos = 0;
		for(int j = 0; j < nface; ++j)  //Now read connectivities, which are triplets of triangle indices
		{
			getline(infile,sline);
			mystringstream.clear();
			mystringstream.str(sline);
			mystringstream >> f2v[pos] >> f2v[pos+1] >> f2v[pos+2];
			pos += 3;
		}
		
		for(int k = 0; k < 3*nface; ++k)  //adjust for C++ 0-based indexing
		{
			f2v[k] -=1;
		}
		
		infile.close();
    
    
    //Write values to mesh's memory arrays
    verts.resize(nvert);
    setCoords(x);

    faces.resize(nface);
    setConnectivities(f2v);

    // clean up
    delete [] x;
    delete [] f2v;

}


/* Write HDF5 mesh */
void Mesh::writeHDF5(const char *fn)
{
    hsize_t dims[2];
    int nvert = numVerts();
    int nface = numFaces();

    double (*x)[3] = new double[nvert][3];
    int (*f2v)[3] = new int[nface][3];

    getCoords(x);
    getConnectivities(f2v);

    hid_t fid = H5Fcreate(fn, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // coordinates
    dims[0] = nvert;
    dims[1] = 3;
    H5LTmake_dataset_double(fid, "X", 2, dims, *x);

    // connectivities
    dims[0] = nface;
    dims[1] = 3;
    H5LTmake_dataset_int(fid, "F2V", 2, dims, *f2v);

    H5Fclose(fid);

    delete [] x;
    delete [] f2v;
}


/* Write the mesh to a Tecplot file in finite-element format
 * Arguments:
 *  fn -- file name */
void Mesh::writeTecplot(const char *fn) 
{
    FILE *file = fopen(fn, "w");

    fprintf(file, "VARIABLES = X, Y, Z\n");
    fprintf(file, "ZONE N=%d E=%d F=FEPOINT ET=TRIANGLE\n", numVerts(), numFaces());

    for (int ivert = 0; ivert < verts.size(); ivert++) {
        Point &vert = verts[ivert];
	fprintf(file, "%15.5E %15.5E %15.5E\n", vert.x[0], vert.x[1], vert.x[2]);
    }

    for (int iface = 0; iface < faces.size(); iface++) { 
        Tri &face = faces[iface];
	fprintf(file, "%9d %9d %9d\n", 
		face.ivert[0] + 1, face.ivert[1] + 1, face.ivert[2] + 1);
    }

    fclose(file);
}


/* Build edge list
 * Arguments:
 *  nvert -- number of vertices
 *  nface -- number of faces
 *  f2v -- face to vertex list
 *  nedge -- number of edges
 *  e2v -- edge to vertex list
 *  e2f -- edge to face list
 * Note:
 *  1. e2v and e2f must be reserved with enough space (3*nface for the 
 *  	first dimension)
 *  2. Each edge has the same orientation as its first neighbor face */
void Mesh_edgeList(int nvert, int nface, const int (*f2v)[3], 
	int &nedge, int (*e2v)[2], int (*e2f)[2])
{
    int NEDGE_MAX = 3*nface;
    int *firstEdge = new int[nvert]; // the first edge connecting to the i-th vertex
    int *nextEdge = new int[NEDGE_MAX];

    // Init
    std::fill(firstEdge, firstEdge + nvert, -1);
    std::fill(nextEdge, nextEdge + NEDGE_MAX, -1);

    // Go over every face
    nedge = 0;
    for (int iface = 0; iface < nface; ++iface) {
	for (int l = 0; l < 3; ++l) {
	    int ivert0 = f2v[iface][l];
	    int ivert1 = f2v[iface][(l+1)%3];

	    if (ivert0 > ivert1) swap(ivert0, ivert1);

	    int p = firstEdge[ivert0];
	    while (p >= 0) {
	        if (e2v[p][0] == ivert0 && e2v[p][1] == ivert1) break;
	        p = nextEdge[p];
	    }

	    if (p < 0) {
		// Create a new edge
		p = nedge;

		e2v[p][0] = ivert0;
		e2v[p][1] = ivert1;

		e2f[p][0] = e2f[p][1] = -1;

		nextEdge[p] = firstEdge[ivert0];
		firstEdge[ivert0] = p;

		nedge++;
	   }

	   assert(e2f[p][0] < 0 || e2f[p][1] < 0);

	   if (e2f[p][0] < 0)
	       e2f[p][0] = iface;
	   else 
	       e2f[p][1] = iface;
        }
    }

    // Make each edge having the same orientation as its first neighbor
    for (int iedge = 0; iedge < nedge; ++iedge) {
        int ivert0 = e2v[iedge][0];
        int ivert1 = e2v[iedge][1];

        int iface0 = e2f[iedge][0];
	if (iface0 < 0) continue;

        for (int l = 0; l < 3; l++) {
            if (f2v[iface0][(l+1)%3] == ivert0 && f2v[iface0][l] == ivert1) {
                swap(e2v[iedge][0], e2v[iedge][1]);
                break;
            }
        }
    }

    // Dealloc temp arrays
    delete []firstEdge;
    delete []nextEdge;
}


/* Build the 1-ring neighbor list
 * Arguments:
 *  nvert, nface, f2v 
 *  nbrs[i] -- neighbor list of the i-th vertex
 * Note:
 *  1. Neighbors are arranged counter-clockwise */
void Mesh_vertOneRingNbrs(int nvert, int nface, const int (*f2v)[3], vector<vector<int> > &nbr)
{
    // Find all edges surrounding every vertex
    // And then store them in a chain list
    vector<int> firstEdge(nvert, -1);
    vector<EdgeWithLnk> edges;

    edges.reserve(3*nface);

    for (int iface = 0; iface < nface; iface++)
    for (int l = 0; l < 3; l++) {
        int ivert = f2v[iface][l];
	edges.push_back( EdgeWithLnk(f2v[iface][(l+1)%3], f2v[iface][(l+2)%3], firstEdge[ivert]) );
	firstEdge[ivert] = edges.size() - 1;
    }

    // Now build neighbor list
    nbr.resize(nvert);

    for (int ivert = 0; ivert < nvert; ivert++) {
	vector<int> &nbr_tmp = nbr[ivert];
	nbr_tmp.clear();
	nbr_tmp.push_back(edges[firstEdge[ivert]].v0());

	while (1) {
	    bool found = false;

	    // Search in both directions 
	    for (int ptr = firstEdge[ivert]; ptr != -1; ptr = edges[ptr].next() ) {

	        int v0 = edges[ptr].v0();
	        int v1 = edges[ptr].v1();

		if (v0 == nbr_tmp.back()) {
		    nbr_tmp.push_back(v1);

		    edges[ptr]._v0 = edges[ptr]._v1 = -1;
		    found = true;
		    break;
		} 
		else if (v1 == nbr_tmp.front()) {
		    nbr_tmp.insert(nbr_tmp.begin(), v0);

		    edges[ptr]._v0 = edges[ptr]._v1 = -1;
		    found = true;
		    break;
		}
	    }

	    if (! found) break;
	}

	// Check that all edges are counted
	for (int ptr = firstEdge[ivert]; ptr != -1; ptr = edges[ptr].next()) {
	    assert( edges[ptr].v0() == -1 && edges[ptr].v1() == -1 );
	}

	// Remove the duplicate front and end points when the neighbors form a closed loop
	if (nbr_tmp.front() == nbr_tmp.back()) nbr_tmp.pop_back();
    }
}


/* Face-to-edge and face-to-face neighbor-list 
 * Arguments:
 *  nvert --
 *  nface --
 *  f2v -- face-to-vert
 *  e2v -- edge-to-vert 
 *  nedge --
 *  e2v -- edge-to-vert
 *  e2f -- edge-to-face 
 *  f2e -- face-to-edge, output
 *  f2f -- face-to-face, output
 * Note:
 *  -- The edge f2e[iface][l] is opposite to f2v[iface][l] 
 *  -- The face f2f[iface][l] is also opposite to f2v[iface][l] */
void Mesh_faceNbrList(int nvert, int nface, const int (*f2v)[3],
		int nedge, const int (*e2v)[2], const int (*e2f)[2], 
		int (*f2e)[3], int (*f2f)[3])
{
    // Init
    for (int iface = 0; iface < nface; iface++) {
        f2e[iface][0] = f2e[iface][1] = f2e[iface][2] = -1;
        f2f[iface][0] = f2f[iface][1] = f2f[iface][2] = -1;
    }

    for (int iedge = 0; iedge < nedge; iedge++) {
	for (int l = 0; l < 2; l++) {
	    int iface = e2f[iedge][l];
	    if (iface < 0) continue;

	    for (int m = 0; m < 3; m++) {
	        int ivert1 = f2v[iface][(m+1)%3];
		int ivert2 = f2v[iface][(m+2)%3];

		if ( (ivert1 == e2v[iedge][0] && ivert2 == e2v[iedge][1]) ||
		     (ivert1 == e2v[iedge][1] && ivert2 == e2v[iedge][0]) ) {
		    f2e[iface][m] = iedge;
		    f2f[iface][m] = e2f[iedge][(l+1)%2];
		    break;
		}
	    } // m
	} // l
    } // iedge
}


/* Reorder face-to-vertex list to make it self-consistent
 * Arguments:
 *  nvert -- number of vertices
 *  nface -- number of faces
 *  f2v -- face-to-vertex list
 * Note:
 *   1. The mesh should be single-connected
 *   2. After orientation, all faces are consistently oriented.  The
 *      orientations are all outward or inward.  */
void Mesh_reorderF2v(int nvert, int nface, int (*f2v)[3])
{
    const int nedgeMax = 3*nface;
    int nedge;
    int (*e2v)[2] = new int[nedgeMax][2];
    int (*e2f)[2] = new int[nedgeMax][2];
    int (*f2e)[3] = new int[nface][3];
    int (*f2f)[3] = new int[nface][3];

    Mesh_edgeList(nvert, nface, f2v, nedge, e2v, e2f);
    Mesh_faceNbrList(nvert, nface, f2v, nedge, e2v, e2f, f2e, f2f);

    int reverse_cnt = 0;    // number of reverted faces

    bool *oriented = new bool[nface];
    std::fill(oriented, oriented + nface, false);

    while (1) {
	// Find the first unoriented face 
        int firstFace = 0;
	while (firstFace < nface && (oriented[firstFace])) firstFace++;
	if (firstFace >= nface) break;    // all faces are handled

        // Flood scheme with firstFace as seed
        vector<int> frontier;
	frontier.push_back(firstFace);
	oriented[firstFace] = true;

	while (!frontier.empty()) {
	    int iface = frontier.back();
	    frontier.pop_back();

	    for (int l = 0; l < 3; ++l) {
		int iface1 = f2f[iface][l]; 
		if (iface1 < 0) continue;
		if (oriented[iface1]) continue;

		bool consistent = false;
		for (int p = 0; p < 3; p++)
		for (int q = 0; q < 3; q++) {
		    if (f2v[iface][p] == f2v[iface1][(q+1)%3] && 
			f2v[iface][(p+1)%3] == f2v[iface1][q]) {
			consistent = true;
			break;
		    }
		} 

		if (! consistent) {
		    swap(f2v[iface1][1], f2v[iface1][2]);
		    reverse_cnt++;
		}
		oriented[iface1] = true;
		frontier.push_back(iface1);
	    }  // l
	}  // while
    } // while

    printf("Re-order face-to-vertex list\n");
    printf("  %d faces have orientation reverted\n", reverse_cnt);

    // Deallocate working arrays
    delete [] e2v;
    delete [] e2f;
    delete [] f2e;
    delete [] f2f;
    delete [] oriented;
}


/* Slice a mesh and calculate the volume in each z-interval
 * Arguments:
 *   zmin, zmax -- min and max z coordinates
 *   nz -- number of intervals
 *   vol_slice -- volume in every interval */
void Mesh_slice(const Mesh &mesh, double zmin, double zmax, int nz, MArray<double,1> &vol_slice)
{
    // Generate cutting planes
    double dz = (zmax - zmin)/nz;
    vector<double> zcut(nz-1);
    for (int iz = 0; iz < nz-1; iz++) zcut[iz] = zmin + (iz+1)*dz;

    // Init
    vol_slice = 0.0;

    for (int ifa = 0; ifa < mesh.numFaces(); ifa++) {
	const Tri &face = mesh.faces[ifa];

	double xtri[3][3];
	for (int l = 0; l < 3; l++) {
	    int ivert = face.ivert[l];
	    m_dcopy(3, mesh.verts[ivert].x, xtri[l]);
	}

	const double znormal[3] = { 0.0, 0.0, 1.0 };
	int nsub;
	int *npsub = new int[nz];
	double (*xsub)[5][3] = new double[nz][5][3];
	slice_triangle(xtri, znormal, zcut, nsub, npsub, xsub);

	for (int isub = 0; isub < nsub; isub++) {
	    double areaSub, xcSub[3];
	    calcPolygonAreaCenter(npsub[isub], xsub[isub], areaSub, xcSub);

	    // Add volume fraction and osmotic in the slab
	    int iz = (int)floor((xcSub[2] - zmin)/dz);
	    if (iz < 0 || iz >= nz) {
	        continue;
	    }

	    // Contribution to the mesh volume
	    double dV = areaSub*(xcSub[1] - mesh.center[1])*face.normal[1];
	    vol_slice(iz) += dV;
	} 

	delete [] npsub;
	delete [] xsub;
    } 

    // Error check
    double vsum = 0.0;
    for (int i = 0; i < nz; i++) vsum += vol_slice(i);

    if (fabs(vsum - mesh.vol) > 1.E-10*mesh.vol) {
	printf("Mesh_slice: mesh is not sliced properly\n");
	printf("    vslice sum = %15.5E\n", vsum);
	printf("    vol        = %15.5E\n", mesh.vol);
    }
}
