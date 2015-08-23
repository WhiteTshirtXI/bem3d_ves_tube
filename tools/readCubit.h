#ifndef _READ_CUBIT
#define _READ_CUBIT

#include "mesh.h"
#include "stokesys.h"

void readCubitMesh(const char *fn, Mesh &mesh);

void readCubitRestart(const char *fn, StokeSys &);

#endif
