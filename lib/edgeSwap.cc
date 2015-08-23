#include "cxxheaders.h"
#include "mathfunc.h"
#include "mesh.h"

/* minAngleOfTri
 * minAngleOfEdge
 * pvalueOfEdge
 * Mesh_edgeSwap */

const double minAngleThresh = M_PI/8.0;

/* The minimum angle of a triangle 
 * Arguments:
 *   x0, x1, x2 -- the coordinates of the three vertices */
double minAngleOfTri(const double *x0, const double *x1, const double *x2)
{
    double xx[3], l0, l1, l2;

    FOR_I3 xx[i] = x1[i] - x0[i];
    l0 = normVec3D(xx);

    FOR_I3 xx[i] = x2[i] - x1[i];
    l1 = normVec3D(xx);

    FOR_I3 xx[i] = x0[i] - x2[i];
    l2 = normVec3D(xx);

    // Make l2 the shortest edge
    if (l2 > l0) swap(l2, l0);
    if (l2 > l1) swap(l2, l1);
    return acos( (l0*l0 + l1*l1 - l2*l2)/(2*l0*l1) );
}


double minAngleOfTri(const Tri &T)
{
    return minAngleOfTri(T.vert[0]->x, T.vert[1]->x, T.vert[2]->x);
}


double minAngleOfEdge(const Edge &edge) 
{
    Tri &T0 = *edge.face[0];
    Tri &T1 = *edge.face[1];

    double th0 = minAngleOfTri(T0);
    double th1 = minAngleOfTri(T1);
    return min(th0, th1);
}


/* Priority value of edge swap, defined as the thMinNew/thMin
 * where thMinNew and thMin are the edge minimal angles after
 * and before the swapping 
 *
 * Before swapping
 *             1
 *           / | \
 *          2  |  3    
 *           \ | /
 *             0
 *
 * After swapping
 *             1
 *           /   \
 *          2-----3    
 *           \   /
 *             0
 * */
double pvalueOfEdge(const Edge &edge)
{
    Point *v0 = edge.vert[0];
    Point *v1 = edge.vert[1];

    Tri *T0 = edge.face[0];
    Tri *T1 = edge.face[1];

    Point *v2;
    for (int l = 0; l < 3; l++) {
        v2 = T0->vert[l];
        if ( v2 != v0 && v2 != v1 ) break;
    }

    Point *v3;
    for (int l = 0; l < 3; l++) {
        v3 = T1->vert[l];
        if ( v3 != v0 && v3 != v1 ) break;
    }

    // Firstly, check that the edge swapping will not flip
    // the normal direction
    double nrml0[3], nrml1[3], nrml0_new[3], nrml1_new[3];
    double detj;
    tri_normal(v0->x, v1->x, v2->x, nrml0, detj);
    tri_normal(v0->x, v3->x, v1->x, nrml1, detj);
    tri_normal(v2->x, v3->x, v1->x, nrml0_new, detj);
    tri_normal(v2->x, v0->x, v3->x, nrml1_new, detj);

    if ( m_ddot(3, nrml0_new, nrml0) <= 0 ||
         m_ddot(3, nrml0_new, nrml1) <= 0 ||
         m_ddot(3, nrml1_new, nrml0) <= 0 ||
         m_ddot(3, nrml1_new, nrml1) <= 0 ) return 0;

    double th = min( minAngleOfTri(v0->x, v1->x, v2->x),
                     minAngleOfTri(v0->x, v1->x, v3->x) );
    double thNew = min( minAngleOfTri(v0->x, v3->x, v2->x),
                        minAngleOfTri(v1->x, v2->x, v3->x) );
    return thNew/th;
}


/* Perform edge swap on a mesh to optimize the mesh quality
 * Arguments:
 *   mesh -- */
void Mesh_edgeSwap(Mesh &mesh)
{
    // Build connectivity information
    int nvert = mesh.numVerts();
    int nface = mesh.numFaces();
    int nedgeMax = nface*3;
    int nedge;

    int (*f2v)[3] = new int[nface][3];
    int (*f2e)[3] = new int[nface][3];
    int (*f2f)[3] = new int[nface][3];
    int (*e2v)[2] = new int[nedgeMax][2];
    int (*e2f)[2] = new int[nedgeMax][2];

    mesh.getConnectivities(f2v);
    Mesh_edgeList(nvert, nface, f2v, nedge, e2v, e2f);
    Mesh_faceNbrList(nvert, nface, f2v, nedge, e2v, e2f, f2e, f2f);

    double *minAngle = new double[nedge];
    double *pval = new double[nedge];

    for (int i = 0; i < nedge; i++) {
        Edge edge;

	edge.vert[0] = &mesh.verts[e2v[i][0]];
	edge.vert[1] = &mesh.verts[e2v[i][1]];

	edge.face[0] = &mesh.faces[e2f[i][0]];
	edge.face[1] = &mesh.faces[e2f[i][1]];

	minAngle[i] = minAngleOfEdge(edge);
	pval[i] = pvalueOfEdge(edge);
    }


    // Do edge swapping iteratively until it can not be done
    for (int iter = 0; iter < nedge; iter++) {
	// Find the edge with the max pvalue
	int ie0 = -1;
	for (int i = 0; i < nedge; i++) {
	    if (minAngle[i] > minAngleThresh) continue;
	    if (pval[i] < 1.0) continue;

	    if (ie0 < 0)
	        ie0 = i;
	    else if (pval[ie0] < pval[i]) 
	    	ie0 = i;
	}

	if (ie0 < 0) break;

	// debug
	printf("%4dth edge has max pvalue\n", ie0);
	// end debug

	// Swap the edge
	// Identify all the vertices, faces, and edges
	int if0 = -1, if1 = -1;
	int iv0 = -1, iv1 = -1, iv2 = -1, iv3 = -1;
	int ie1 = -1, ie2 = -1, ie3 = -1, ie4 = -1;
	
	if0 = e2f[ie0][0];
	if1 = e2f[ie0][1];

	iv0 = e2v[ie0][0];
	iv1 = e2v[ie0][1];
	iv2 = f2v[if0][0] + f2v[if0][1] + f2v[if0][2] - iv0 - iv1;
	iv3 = f2v[if1][0] + f2v[if1][1] + f2v[if1][2] - iv0 - iv1;

	for (int l = 0; l < 3; l++) {
	    int ie = f2e[if0][l];
	    if ( ( e2v[ie][0] == iv2 && e2v[ie][1] == iv0 ) ||
	         ( e2v[ie][0] == iv0 && e2v[ie][1] == iv2 ) ) {
	        ie1 = ie; 
	    } 

	    if ( ( e2v[ie][0] == iv1 && e2v[ie][1] == iv2 ) ||
                 ( e2v[ie][0] == iv2 && e2v[ie][1] == iv1 ) ) {
		ie2 = ie;
	    }

	    ie = f2e[if1][l];
	    if ( ( e2v[ie][0] == iv0 && e2v[ie][1] == iv3 ) ||
	         ( e2v[ie][0] == iv3 && e2v[ie][1] == iv0 ) ) {
	        ie3 = ie; 
	    } 

	    if ( ( e2v[ie][0] == iv1 && e2v[ie][1] == iv3 ) ||
	         ( e2v[ie][0] == iv3 && e2v[ie][1] == iv1 ) ) {
	        ie4 = ie; 
	    } 
	}

	// Check consistency
	assert( if0 >= 0 && if1 >= 0 &&
		iv0 >= 0 && iv1 >= 0 && iv2 >= 0 && iv3 >= 0 &&
		ie0 >= 0 && ie1 >= 0 && ie2 >= 0 && ie3 >= 0 && ie4 >= 0);

	// Update e2v
	e2v[ie0][0] = iv2;
	e2v[ie0][1] = iv3;
	
	// Update f2v and f2e
	f2v[if0][0] = iv2;
	f2v[if0][1] = iv3;
	f2v[if0][2] = iv1;

	f2e[if0][0] = ie4;
	f2e[if0][1] = ie2;
	f2e[if0][2] = ie0;

	f2v[if1][0] = iv2;
	f2v[if1][1] = iv0;
	f2v[if1][2] = iv3;

	f2e[if1][0] = ie3;
	f2e[if1][1] = ie0;
	f2e[if1][2] = ie1;

	// Update e2f
	for (int l = 0; l < 2; l++) {
	    if (e2f[ie1][l] == if0) e2f[ie1][l] = if1;
	    if (e2f[ie4][l] == if1) e2f[ie4][l] = if0;
	}

	// Update mesh
	Tri *face = &mesh.faces[if0];
	FOR_I3 {
	    face->ivert[i] = f2v[if0][i];
	    face->vert[i] = &mesh.verts[face->ivert[i]];
	}

	face = &mesh.faces[if1];
	FOR_I3 {
	    face->ivert[i] = f2v[if1][i];
	    face->vert[i] = &mesh.verts[face->ivert[i]];
	}

	// Update edge minimum angle and pvalue
	int edge_update[5] = { ie0, ie1, ie2, ie3, ie4 };
	for (int p = 0; p < 5; p++) {
	    int i = edge_update[p];

	    Edge edge;

	    edge.vert[0] = &mesh.verts[e2v[i][0]];
	    edge.vert[1] = &mesh.verts[e2v[i][1]];

	    edge.face[0] = &mesh.faces[e2f[i][0]];
	    edge.face[1] = &mesh.faces[e2f[i][1]];

	    minAngle[i] = minAngleOfEdge(edge);
	    pval[i] = pvalueOfEdge(edge);
	}
    } // while

    delete [] f2v;
    delete [] f2e;
    delete [] f2f;

    delete [] e2v;
    delete [] e2f;

    delete [] minAngle;
    delete [] pval;
}
