#include "cxxheaders.h"
#include "readCubit.h"

int main(int argc, char **agrv)
{
    StokeSys stks;
    char fn[256];

    readCubitRestart(agrv[1], stks);

    printf("New restart file in HDF5 format =? ");
    scanf("%s", fn);
    stks.writeRestart(fn);
}
