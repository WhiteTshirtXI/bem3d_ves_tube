#include "cxxheaders.h"
#include "mesh.h"
#include "cell.h"
#include "rigid.h"
#include "wall.h"
#include "stokesys.h"
#include "ewald.h"
#include "netcdf.h"
#include "readCubit.h"

/* readCubitMesh
 * readCubitRestart */

/* Read cubit mesh
 * Arguments:
 *  fn -- Cubit file name 
 *  mesh -- the mesh */
void readCubitMesh(const char *fn, Mesh &mesh)
{
    int ncid, varid, dimids[2];
    size_t dimlens[2];
    char var_name[512];

    nc_open(fn, NC_NOWRITE, &ncid);

    // Read coordinates
    MArray<double,2> x;

    if (nc_inq_varid(ncid, "coord", &varid) == NC_NOERR) {
    	nc_inq_vardimid(ncid, varid, dimids);
	nc_inq_dimlen(ncid, dimids[0], &dimlens[0]);
	nc_inq_dimlen(ncid, dimids[1], &dimlens[1]);

	MArray<double,2> xbuf(dimlens[0], dimlens[1]);
	nc_get_var_double(ncid, varid, xbuf.data());

	x.resize(dimlens[1], 3);
	x = 0.0;
	for (int i = 0; i < dimlens[1]; i++)
	for (int j = 0; j < dimlens[0]; j++) {
	    x(i,j) = xbuf(j,i);
	}
    }
    else if (nc_inq_varid(ncid, "coordx", &varid) == NC_NOERR) {
    	nc_inq_vardimid(ncid, varid, dimids);
	nc_inq_dimlen(ncid, dimids[0], &dimlens[0]);

	MArray<double,1> xbuf(dimlens[0]);
	xbuf = 0.0;
	nc_get_var_double(ncid, varid, xbuf.data());

	MArray<double,1> ybuf(dimlens[0]);
	ybuf = 0.0;
	if (nc_inq_varid(ncid, "coordy", &varid) == NC_NOERR) {
	    nc_get_var_double(ncid, varid, ybuf.data());
	}

	MArray<double,1> zbuf(dimlens[0]);
	zbuf = 0.0;
	if (nc_inq_varid(ncid, "coordz", &varid) == NC_NOERR) {
	    nc_get_var_double(ncid, varid, zbuf.data());
	}

	x.resize(dimlens[0], 3);
	for (int i = 0; i < dimlens[0]; i++) {
	    x(i,0) = xbuf(i);
	    x(i,1) = ybuf(i);
	    x(i,2) = zbuf(i);
	}
    }


    // Read connectivities
    vector<int> f2v;

    for (int p = 1; ; p++) {
        sprintf(var_name, "connect%d", p);
	if (nc_inq_varid(ncid, var_name, &varid) != NC_NOERR) break;

	nc_inq_vardimid(ncid, varid, dimids);
	nc_inq_dimlen(ncid, dimids[0], &dimlens[0]);
	nc_inq_dimlen(ncid, dimids[1], &dimlens[1]);
	assert(dimlens[1] == 3);

	MArray<int,1> f2v_tmp(dimlens[0]*dimlens[1]);
        nc_get_var_int(ncid, varid, f2v_tmp.data());

	f2v.insert(f2v.end(), &f2v_tmp(0), &f2v_tmp(0) + dimlens[0]*dimlens[1]);
    }

    // Copy contents to mesh
    int nvert = x.size(0);
    int nface = f2v.size()/3;

    mesh.verts.resize(nvert);
    mesh.faces.resize(nface);

    for (int ivert = 0; ivert < nvert; ivert++) {
        Point &vert = mesh.verts[ivert];
        FOR_J3 vert.x[j] = x(ivert,j);
    }

    for (int iface = 0; iface < nface; iface++) {
        Tri &face = mesh.faces[iface];
        FOR_J3 face.ivert[j] = f2v[3*iface+j] - 1;
    }
}


/*
void readCubitRestart(const char *fn, StokeSys &stks)
{
    int ncid, varid, dimid;
    size_t dimlen;
    int ncell, nrigid, nwall;
    char token[256];
    int status;

    printf("Read restart file %s\n", fn);

    nc_open(fn, NC_NOWRITE, &ncid);

    nc_inq_varid(ncid, "EWALD_LX", &varid);
    nc_get_var_double(ncid, varid, &ewald::L[0]);

    nc_inq_varid(ncid, "EWALD_LY", &varid);
    nc_get_var_double(ncid, varid, &ewald::L[1]);

    nc_inq_varid(ncid, "EWALD_LZ", &varid);
    nc_get_var_double(ncid, varid, &ewald::L[2]);

    nc_inq_varid(ncid, "TIME", &varid);
    nc_get_var_double(ncid, varid, &stks.time);

    nc_inq_varid(ncid, "LT", &varid);
    nc_get_var_int(ncid, varid, &stks.lt);

    nc_inq_varid(ncid, "VBKG_X", &varid);
    nc_get_var_double(ncid, varid, &stks.vbkg[0]);
    
    nc_inq_varid(ncid, "VBKG_Y", &varid);
    nc_get_var_double(ncid, varid, &stks.vbkg[1]);

    nc_inq_varid(ncid, "VBKG_Z", &varid);
    nc_get_var_double(ncid, varid, &stks.vbkg[2]);

    printf("  L = %.3f  %.3f  %.3f\n", ewald::L[0], ewald::L[1], ewald::L[2]);
    printf("  time = %.3f\n", stks.time);
    printf("  lt = %d\n", stks.lt);
    printf("  vbkg = %.3f  %.3f  %.3f\n", stks.vbkg[0], stks.vbkg[1], stks.vbkg[2]);

    ncell = 0;
    status = nc_inq_dimid(ncid, "NCELL", &dimid);
    if (status == NC_NOERR) {
	nc_inq_dimlen(ncid, dimid, &dimlen);
	ncell = dimlen;
    } 

    nrigid = 0;
    status = nc_inq_dimid(ncid, "NRIGID", &dimid);
    if (status == NC_NOERR) {
	nc_inq_dimlen(ncid, dimid, &dimlen);
	nrigid = dimlen;
    } 

    nwall = 0;
    status = nc_inq_dimid(ncid, "NWALL", &dimid);
    if (status == NC_NOERR) {
	nc_inq_dimlen(ncid, dimid, &dimlen);
	nwall = dimlen;
    } else

    printf("  ncell = %d\n", ncell);
    printf("  nrigid = %d\n", nrigid);
    printf("  nwall = %d\n", nwall);

    stks.cells.resize(ncell);
    stks.rigids.resize(nrigid);
    stks.walls.resize(nwall);

    // Read cells
    for (int icell = 0; icell < ncell; icell++) {
        Cell &cell = stks.cells[icell];
        int nvert, nface;

	sprintf(token, "CELL_%d_NVERT", icell);
	nc_inq_dimid(ncid, token, &dimid);
	nc_inq_dimlen(ncid, dimid, &dimlen);
	nvert = dimlen;

	sprintf(token, "CELL_%d_NFACE", icell);
	nc_inq_dimid(ncid, token, &dimid);
	nc_inq_dimlen(ncid, dimid, &dimlen);
	nface = dimlen;

        double *x = new double[3*nvert];
	int *f2v = new int[3*nface];

	sprintf(token, "CELL_%d_X", icell);
	nc_inq_varid(ncid, token, &varid);
	nc_get_var_double(ncid, varid, x);

	sprintf(token, "CELL_%d_F2V", icell);
	nc_inq_varid(ncid, token, &varid);
	nc_get_var_int(ncid, varid, f2v);

	cell.verts.resize(nvert);
	cell.setCoords(x);

	cell.faces.resize(nface);
	cell.setConnectivities(f2v);

	sprintf(token, "CELL_%d_XREF", icell);
	status = nc_inq_varid(ncid, token, &varid);

	if (status == NC_NOERR) {
	    nc_get_var_double(ncid, varid, x);

	    cell.cellRef = new Cell();

	    cell.cellRef->verts.resize(nvert);
	    cell.cellRef->setCoords(x);

	    cell.cellRef->faces.resize(nface);
	    cell.cellRef->setConnectivities(f2v);
	}

	// dealloc temp arrays
	delete [] x;
	delete [] f2v;
    } 


    // Read rigids
    for (int irigid = 0; irigid < nrigid; irigid++) {
        Rigid &rigid = stks.rigids[irigid];
        int nvert, nface;

	sprintf(token, "RIGID_%d_NVERT", irigid);
	nc_inq_dimid(ncid, token, &varid);
	nc_inq_dimlen(ncid, varid, &dimlen);
	nvert = dimlen;

	sprintf(token, "RIGID_%d_NFACE", irigid);
	nc_inq_dimid(ncid, token, &varid);
	nc_inq_dimlen(ncid, varid, &dimlen);
	nface = dimlen;

        double *x = new double[3*nvert];
	int *f2v = new int[3*nface];

	sprintf(token, "RIGID_%d_X", irigid);
	nc_inq_varid(ncid, token, &varid);
	nc_get_var_double(ncid, varid, x);

	sprintf(token, "RIGID_%d_F2V", irigid);
	nc_inq_varid(ncid, token, &varid);
	nc_get_var_int(ncid, varid, f2v);

	rigid.verts.resize(nvert);
	rigid.setCoords(x);

	rigid.faces.resize(nface);
	rigid.setConnectivities(f2v);

	// dealloc temp arrays
	delete [] x;
	delete [] f2v;
    }


    // Read walls
    for (int iwall = 0; iwall < nwall; ++iwall) {
        Wall &wall = stks.walls[iwall];
        int nvert, nface;

	sprintf(token, "WALL_%d_NVERT", iwall);
	nc_inq_dimid(ncid, token, &varid);
	nc_inq_dimlen(ncid, varid, &dimlen);
	nvert = dimlen;

	sprintf(token, "WALL_%d_NFACE", iwall);
	nc_inq_dimid(ncid, token, &varid);
	nc_inq_dimlen(ncid, varid, &dimlen);
	nface = dimlen;

        double *x = new double[3*nvert];
	int *f2v = new int[3*nface];

	sprintf(token, "WALL_%d_X", iwall);
	nc_inq_varid(ncid, token, &varid);
	nc_get_var_double(ncid, varid, x);

	sprintf(token, "WALL_%d_F2V", iwall);
	nc_inq_varid(ncid, token, &varid);
	nc_get_var_int(ncid, varid, f2v);

	wall.verts.resize(nvert);
	wall.setCoords(x);

	wall.faces.resize(nface);
	wall.setConnectivities(f2v);

	// surface force densities
	wall.f = new double[nvert][3];
	sprintf(token, "WALL_%d_F", iwall);
	status = nc_inq_varid(ncid, token, &varid);

	if (status == NC_NOERR) nc_get_var_double(ncid, varid, *wall.f);

	// dealloc temp arrays
	delete []x;
	delete []f2v;
    }  // iwall

    nc_close(ncid);
} */
