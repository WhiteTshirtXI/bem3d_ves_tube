include ../Makefile.in

LIB = $(COMMON_LIB) $(HDF5_LIB) $(FFTW2_LIB) $(PETSC_LIB)

BINS = $(patsubst %.cc,%,$(wildcard *.cc))

EXES =
#EXES += slipVel
###EXES += centerPlate
###EXES += velField
###EXES += cellDens
###EXES += pressGrad
#EXES += tracer_shot
#EXES += rename_fields
##EXES += rigidTraj.02
#EXES += analyze_probe
#EXES += analyze_wall
###EXES += tracer_penetrate
#EXES += analyze_rigid
#EXES += resetTime
#EXES += add_probe_tracer
#EXES += add_rigid
#EXES += channel
#EXES += analyze_tracer
#EXES += analyze_cell
EXES += vesinitcond
EXES += vchannel
#EXES += fahraeus
#EXES += addctc
#EXES += addctc_findmode
#EXES += mandeletecell
#EXES += rigcellswap
#EXES += ctcCOM
#EXES += ctcCOMR
#EXES += wallclosest
#EXES += add_rigid_spann
#EXES+=test
#EXES += addrspann

all : $(EXES)

lib :
	make -C $(WORK_DIR)/lib

$(EXES) : % : %.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIB)

clean :
	-rm -f $(EXES) $(BINS) *.o *.a core

# Dependency
depend .depend :
	$(CXX) $(INCLUDE) -MM $(wildcard *.cc) > .depend

include .depend
