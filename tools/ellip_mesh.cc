#include "cxxheaders.h"
#include "mesh.h"
#include "mathfunc.h"
#include "miscUtils.h"

void calcRelaxVel(Mesh &, MArray<double,2> &);

int main(int argc, char **argv)
{
    Mesh mesh;

    // Squeeze the sphere
    mesh.makeSphere(6);
    for (int ivert = 0; ivert < mesh.numVerts(); ivert++) {
        Point &vert = mesh.verts[ivert];
	vert.x[2] *= 0.3333;
    }

    mesh.setInternalPointers();

/*
    // Relaxation
    int nvert = mesh.numVerts();
    int nface = mesh.numFaces();

    for (int iter = 0; iter < 100; iter++) {
        MArray<double,2> v;

	mesh.updateGeometry();
	calcRelaxVel(mesh, v);

	for (int ivert = 0; ivert < nvert; ivert++) {
	    Point &vert = mesh.verts[ivert];
	    FOR_D3 vert.x[d] += 0.1*v(ivert,d);
	}
    }
*/

    miscUtils::mkdir("D");
    mesh.writeTecplot("D/ellipsoid.dat");
    mesh.writeHDF5("D/ellipsoid.h5");
}


/*
void calcRelaxVel(Mesh &mesh, MArray<double,2> &v)
{
    int nvert = mesh.numVerts();
    int nface = mesh.numFaces();
    int nedge = mesh.numEdges();

    // Init
    v.resize(nvert,3);
    v = 0.0;

    // Vertex valence
    MArray<int,1> valence(nvert);
    valence = 0;
    for (int iface = 0; iface < nface; iface++) {
        for (int l = 0; l < 3; l++) {
            int ivert = mesh.faces[iface].ivert[l];
            valence(ivert)++;
        }
    }

    // Calculate sprint constant
    MArray<double,1> sk(nvert);	

    sk = 1.0;	// Uniform sprint constant 

    // Weighted by the local curvature and 1/(local average triangle area)
    for (int ivert = 0; ivert < nvert; ivert++) {
	sk(ivert) *= pow(fabs(mesh.vertH(ivert)), 0.5);
	sk(ivert) *= pow(mesh.vertArea(ivert)/valence(ivert), 1.0);
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
	Edge &edge = mesh.edges[ied];

	int iv0 = edge.ivert[0];
	int iv1 = edge.ivert[1];

	double xx[3];
	for (int d = 0; d < 3; d++) 
	    xx[d] = mesh.verts[iv1].x[d] - mesh.verts[iv0].x[d];

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
	Point &vert = mesh.verts[ivert];

	double vn = m_ddot(3, &v(ivert,0), &mesh.vertNrml(ivert,0));
	m_daxpy(3, -vn, &mesh.vertNrml(ivert,0), &v(ivert,0));
    } 
} */
