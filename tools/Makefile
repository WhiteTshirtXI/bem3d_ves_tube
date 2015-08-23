include ../Makefile.in

LIB = $(COMMON_LIB) $(FFTW2_LIB) $(PETSC_LIB) $(HDF5_LIB) -lnetcdf

BINS = $(patsubst %.cc,%,$(wildcard *.cc))

EXES = 
#EXES += restart_2_hdf5
EXES += convert_cubit
#EXES += ellip_mesh

all : $(EXES)

lib :
	make -C $(WORK_DIR)/lib

$(EXES) : % : %.o readCubit.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB) 

clean :
	-rm -f $(EXES) $(BINS) *.o core

# Dependency
depend .depend :
	$(CXX) $(INCLUDE) -MM $(wildcard *.cc) > .depend

include .depend
