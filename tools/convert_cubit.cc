// Convert a mesh generated Cubit to HDF5 and Tecplot format
#include "cxxheaders.h"
#include "mesh.h"
#include "readCubit.h"

int main(int argc, char **argv)
{
    char fn[256];
    Mesh mesh;

    strcpy(fn, argv[1]);
//    strcpy(fn, "/home/hong/tmp/cavity_flow.g");
    readCubitMesh(fn, mesh);
    mesh.setInternalPointers();

    // Information
    int nvert = mesh.numVerts();
    int nface = mesh.numFaces();
    printf("Raw mesh: nvert, nface = %d, %d\n", nvert, nface);

    // Reorder f2v list
    {
	int (*f2v)[3] = new int[nface][3];
	mesh.getConnectivities(f2v);
	Mesh_reorderF2v(nvert, nface, f2v);
	mesh.setConnectivities(f2v);
	delete [] f2v;
    }

    // Test if the surface is closed
    bool closed = true;
    for (int iedge = 0; iedge < mesh.numEdges(); iedge++) {
        Edge &edge = mesh.edges[iedge];
	if (edge.iface[0] < 0 || edge.iface[1] < 0) {
	    closed = false;
	    break;
	}
    }

    // Ensure positive volume if the mesh is indeed closed
    if (closed) {
	mesh.setInternalPointers();
	mesh.updateGeometry();

	if (mesh.vol < 0) {
	    printf("Negative volume, revert all faces\n");

	    for (int iface = 0; iface < nface; ++iface) {
		Tri &face = mesh.faces[iface];
		swap(face.ivert[1], face.ivert[2]);
	    }
	}
    } else {
        printf("Surface is not closed\n");
    }

//    // Consolidate duplicate vertices
//    MArray<double,2> x;
//    MArray<int,2> f2v;
//
//    mesh.getCoords(x);
//    mesh.getConnectivities(f2v);
//
//    MArray<int,1> p2p(nvert);
//    for (int i = 0; i < nvert; i++) p2p(i) = i;
//
//    for (int i = 0; i < nvert; i++) {
//        const double eps = 1.E-10;
//
//	for (int j = i-1; j >= 0; j--) {
//	    if ( fabs(x(i,0) - x(j,0)) < eps && 
//	         fabs(x(i,1) - x(j,1)) < eps && 
//		 fabs(x(i,2) - x(j,2)) < eps ) {
//	        p2p(i) = p2p(j);
//		break;
//            }
//        }
//    }
//
//
//    int nvert_new = 0;
//
//    MArray<double,2> x_new(nvert,3);
//    MArray<int,1> id_new(nvert);	// vertex index after compression
//
//    for (int ivert = 0; ivert < nvert; ivert++) {
//        if (p2p(ivert) == ivert) {
//	    id_new(ivert) = nvert_new;
//	    FOR_J3 x_new(nvert_new,j) = x(ivert,j);
//	    nvert_new++;
//        }
//        else
//	    id_new(ivert) = id_new(p2p(ivert));
//    }
//
//    MArray<int,2> f2v_new(nface,3);
//    for (int iface = 0; iface < nface; iface++) {
//        FOR_J3 f2v_new(iface,j) = id_new(f2v(iface,j));
//    }
//
//    if (nvert > nvert_new) {
//        mesh.verts.resize(nvert_new);
//
//	mesh.setCoords(x_new);
//	mesh.setConnectivities(f2v_new);
//
//	printf("Consolidate %d duplicated points\n", nvert - nvert_new);
//    }


    // Post process
    mesh.setInternalPointers();
    mesh.updateGeometry();

    // Output
    // Extract file name head
    char fn_head[256];

    {
        int n = strlen(fn);
	for (; n > 0; n--) {
	    if (fn[n] == '.') break;
        }

	if (n > 0) {
	    for (int i = 0; i < n; i++) {
		fn_head[i] = fn[i];
	    }
	    fn_head[n] = '\0';
        }
	else {
	    strcpy(fn_head, "1");
	}
    }

    sprintf(fn, "%s.h5", fn_head);
    printf("Output HDF5 file   : %s\n", fn);
    mesh.writeHDF5(fn);

    sprintf(fn, "%s.dat", fn_head);
    printf("Output Tecplot file: %s\n", fn);
    mesh.writeTecplot(fn);

}
