#include "cxxheaders.h"
#include "marray.h"
#include "point.h"
#include "tri.h"
#include "mesh.h"
#include "mathfunc.h"

/* void Mesh::makeIcosahedron 
 * void Mesh::makeSphere 
 * void Mesh::makeBiConcave 
 * void Mesh::makeOblate 
 * void Mesh::makeProlate */

// Create an icosahedron of edge length 2
void Mesh::makeIcosahedron()
{
    verts.resize(0);
    faces.resize(0);

    // Create vertices
    double t = 0.5*(1.0 + sqrt(5.0)), s = 1.0;

    verts.push_back(Point(t, s, 0));
    verts.push_back(Point(-t, s, 0));
    verts.push_back(Point(t, -s, 0));
    verts.push_back(Point(-t, -s, 0));

    verts.push_back(Point(s, 0, t));
    verts.push_back(Point(s, 0, -t));
    verts.push_back(Point(-s, 0, t));
    verts.push_back(Point(-s, 0, -t));

    verts.push_back(Point(0, t, s));
    verts.push_back(Point(0, -t, s));
    verts.push_back(Point(0, t, -s));
    verts.push_back(Point(0, -t, -s));

    // Create faces
    for (int i0 = 0; i0 < verts.size(); i0++)
    for (int i1 = i0 + 1; i1 < verts.size(); i1++)
    for (int i2 = i1 + 1; i2 < verts.size(); i2++) {
	double xx[3];
	double r01, r12, r20;

	r01 = distCart3D(verts[i0].x, verts[i1].x);
	r12 = distCart3D(verts[i1].x, verts[i2].x);
	r20 = distCart3D(verts[i2].x, verts[i0].x);

	if (fabs(r01 - 2.0) < 0.01 
	 && fabs(r12 - 2.0) < 0.01 
	 && fabs(r20 - 2.0) < 0.01)
	    faces.push_back( Tri(i0, i1, i2) );
    } // i0, i1, i2

    // There should be 12 vertices and 20 faces
    assert(verts.size() == 12 && faces.size() == 20);

    // Make the face normal pointing outward
    for (vector<Tri>::iterator T = faces.begin(); T < faces.end(); T++) {
	int i0 = T->ivert[0], i1 = T->ivert[1], i2 = T->ivert[2];

	double xx1[3], xx2[3], normal[3], center[3];
	
	FOR_J3 {
	    xx1[j] = verts[i1].x[j] - verts[i0].x[j];
	    xx2[j] = verts[i2].x[j] - verts[i0].x[j];
	    center[j] = 1.0/3.0*(verts[i0].x[j] + verts[i1].x[j] + verts[i2].x[j]);
	}
	cross_product(xx1, xx2, normal);

	if (m_ddot(3, normal, center) < 0) {
	    T->ivert[0] = i1;    
	    T->ivert[1] = i0;    
	    T->ivert[2] = i2;
	}
    }
}


/* Tesselate a sphere by dividing an icosahedron and then
 * project the surface to the unit sphere
 * Arguments:
 *   N -- number of divisions */
void Mesh::makeSphere(const int N)
{
    makeIcosahedron();

    int nfaceOld = faces.size();
    MArray<int,3> f2v(nfaceOld,N+1,N+1);

    // Populate f2v
    for (int iface = 0; iface < nfaceOld; iface++) {
	Tri &face = faces[iface];

	// Add corner points
	f2v(iface,0,0) = face.ivert[0];
	f2v(iface,N,0) = face.ivert[1];
	f2v(iface,0,N) = face.ivert[2];

	// Add edge points
	for (int l = 0; l < 3; l++) {
	    int iBgn = face.ivert[l], iEnd = face.ivert[(l+1)%3];
	    // Do not use reference here, as the memory of verts[iBgn] and verts[iEnd] 
	    // can be changed during the process
	    Point vBgn = verts[iBgn];	
	    Point vEnd = verts[iEnd];

	    if (iBgn < iEnd) {	// Avoid redundancy 
		// Find the triangular element that shares the edge
		int iface1, l1;
		for (iface1 = 0; iface1 < nfaceOld; ++iface1) {
		    Tri &face1 = faces[iface1];
		    for (l1 = 0; l1 < 3; l1++)
			if (face1.ivert[l1] == iEnd && face1.ivert[(l1+1)%3] == iBgn) break;
		    if (l1 < 3) break;
		}
		assert(iface1 < nfaceOld);
		
		for (int m = 1; m < N; ++m) {
		    Point vNew;

		    double s = m*1.0/N;
		    FOR_I3 vNew.x[i] = (1.0 - s)*vBgn.x[i] + s*vEnd.x[i];

		    verts.push_back(vNew);
		    int ivertNew = verts.size() - 1;

		    switch (l) {
			case 0:
			    f2v(iface,m,0) = ivertNew;
			    break;
			case 1:
			    f2v(iface,N-m,m) = ivertNew;
			    break;
			case 2:
			    f2v(iface,0,N-m) = ivertNew;
			    break;
		    }

		    switch(l1) {
			case 0:
			    f2v(iface1,N-m,0) = ivertNew;
			    break;
			case 1:
			    f2v(iface1,m,N-m) = ivertNew;
			    break;
			case 2:
			    f2v(iface1,0,m) = ivertNew;
			    break;
		    }
	        } // m
	    }
	} // l

	// Add interior points
	for (int i = 1; i < N; i++)
	for (int j = 1; j < N-i; j++) {
	    Point vNew;
	    Point &vert0 = verts[face.ivert[0]];
	    Point &vert1 = verts[face.ivert[1]];
	    Point &vert2 = verts[face.ivert[2]];
	    double s1 = double(i)/N;
	    double s2 = double(j)/N;
	    
	    FOR_K3 vNew.x[k] = (1.0 - s1 - s2)*vert0.x[k] + s1*vert1.x[k] + s2*vert2.x[k];
	    verts.push_back(vNew);
	    f2v(iface,i,j) = verts.size() - 1;
	} // i j
    } // iface

    // Generate triangles
    faces.clear();    // delete all old triangles
    for (int iface = 0; iface < nfaceOld; iface++) {
	for (int i = 0; i < N; i++)
	for (int j = 0; j < N-i; j++) {
	    faces.push_back( Tri(f2v(iface,i,j), f2v(iface,i+1,j), f2v(iface,i,j+1)) );
	    if (i + j + 2 <= N)
		faces.push_back( Tri(f2v(iface,i+1,j), f2v(iface,i+1,j+1), f2v(iface,i,j+1)) );
	} // i, j
    } // iface


    // Project points to the unit sphere
    for (int i = 0; i < verts.size(); i++) {
	double rad = normVec3D(verts[i].x);
	m_dscal(3, 1.0/rad, verts[i].x);
    }  // i
}


/* Make a biconcave shaped cell
 * Arguments:
 *   N -- number of subdivisions */
void Mesh::makeBiConcave(int N)
{
    makeSphere(N);

    // Map the sphere into bi-concaved shape
    const double alpha = 1.3858189;

    for (int i = 0; i < verts.size(); i++) {
        double th = acos(verts[i].x[2]);
	double phi = atan2(verts[i].x[1], verts[i].x[0]);

	double costh = cos(th);
	double sinth = sin(th);
	double sinth2 = sinth*sinth;
	double sinth4 = sinth2*sinth2;

	double z = 0.5*alpha*(0.207 + 2.003*sinth2 - 1.123*sinth4)*costh;
	double r = alpha*sinth;

	verts[i].x[0] = r*cos(phi);
	verts[i].x[1] = r*sin(phi);
	verts[i].x[2] = z;
    } // i
}


/* Make an oblate spheroid
 * Arguments:
 *   N -- number of subdivisions
 *   volRat -- targeted reduced volume 
 * Note:
 *   -- surface area = 4*PI */
void Mesh::makeOblate(int N, double volRat)
{
    assert(volRat > 0.0 && volRat < 1.0);

    makeSphere(N);

    int nvert = numVerts();
    MArray<double,2> xref(nvert,3);
    getCoords(xref);

    // L -- the length of the shorter axis
    // Use linear search to find L
    double L1 = 0.0;
    double volRat1 = 0.0;

    double L0 = 1.0;
    double volRat0 = 1.0;

    const int MAX_ITER = 20;
    for (int iter = 0; iter < MAX_ITER; iter++) {
	double L2 = L0 + (L1 - L0)*(volRat - volRat0)/(volRat1 - volRat0);
	double volRat2;

	for (int ivert = 0; ivert < nvert; ivert++) {
	    verts[ivert].x[0] = xref(ivert,0);
	    verts[ivert].x[1] = xref(ivert,1);
	    verts[ivert].x[2] = xref(ivert,2)*L2;
	}

	updateAVC();
	double rad = sqrt(area/(4*M_PI));
	volRat2 = vol/(4.0/3.0*M_PI*rad*rad*rad);

	// Normalization
	for (int ivert = 0; ivert < nvert; ivert++) {
	    m_dscal(3, 1.0/rad, verts[ivert].x);
	}

	// printf("L = %6.3f volRat = %6.3f\n", L2, volRat2);

	if (volRat2 > volRat) {
	    L0 = L2;
	    volRat0 = volRat2;
	}
	else {
	    L1 = L2;
	    volRat1 = volRat2;
	}

	// Convergece test
	if (fabs(volRat2 - volRat) < 0.001) break;
    } // iter
}


/* Make an prolate spheroid
 * Arguments:
 *   N -- number of subdivisions
 *   volTar -- targeted reduced volume 
 * Note:
 *   -- surface area = 4*PI */
void Mesh::makeProlate(int N, double volRat)
{
    assert(volRat > 0.0 && volRat < 1.0);

    makeSphere(N);

    int nvert = numVerts();
    MArray<double,2> xref(nvert,3);
    getCoords(xref);

    // L -- the length of the shorter axis
    // Use linear search to find L
    double L1 = 0.0;
    double volRat1 = 0.0;

    double L0 = 1.0;
    double volRat0 = 1.0;

    const int MAX_ITER = 50;

    for (int iter = 0; iter < MAX_ITER; iter++) {
	double L2 = L0 + (L1 - L0)*(volRat - volRat0)/(volRat1 - volRat0);
	double volRat2;

	for (int ivert = 0; ivert < nvert; ivert++) {
	    verts[ivert].x[0] = xref(ivert,0)*L2;
	    verts[ivert].x[1] = xref(ivert,1)*L2;
	    verts[ivert].x[2] = xref(ivert,2);
	}

	updateAVC();
	double rad = sqrt(area/(4*M_PI));
	volRat2 = vol/(4.0/3.0*M_PI*rad*rad*rad);

	// Normalization
	for (int ivert = 0; ivert < nvert; ivert++) {
	    m_dscal(3, 1.0/rad, verts[ivert].x);
	}

	// printf("L = %6.3f volRat = %6.3f\n", L2, volRat2);

	if (volRat2 > volRat) {
	    L0 = L2;
	    volRat0 = volRat2;
	}
	else {
	    L1 = L2;
	    volRat1 = volRat2;
	}

	if (fabs(volRat2 - volRat) < 0.001) break;
    } // iter
}
