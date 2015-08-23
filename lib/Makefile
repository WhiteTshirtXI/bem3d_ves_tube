include ../Makefile.in

OJBS =
OBJS += point.o tri.o edge.o mesh.o sphere_tess.o edgeSwap.o
OBJS += bspline.o geom_oper.o quadrature.o membrane.o 
OBJS += ewald.o ewald_exact.o ewald_phys.o ewald_four.o collision.o
OBJS += srclist.o nbrlist.o
OBJS += cell.o rigid.o wall.o redcell.o
OBJS += greenfunc.o
#OBJS += bonds.o bondcollection.o
#OBJS +=  stokesys.o stokesys_flowsolver.o \
#	stokesys_cell_contact.o stokesys_collision_force.o
#OBJS +=  shearsys.o shearsys_flowsolver.o \
#	shearsys_cell_contact.o shearsys_collision_force.o
#OBJS +=  rbcsys.o rbcsys_bie.o 
#OBJS +=  vesys.o vesys_bie.o vesys_project.o
OBJS +=  veswall.o veswall_wallBie.o veswall_vesicleBie.o veswall_project.o
OBJS += debugfunc.o mympi.o mypetsc.o
OBJS += param.o miscUtils.o domainDecomp.o
OBJS += subdiv.o vesicle.o 
OBJS += correlation.o



all lib : libBIE.a

libBIE.a : $(OBJS)
	@echo OBJS = $(OBJS)
	$(AR) $@ $?
	$(RANLIB) $@

clean :
	rm -f *.o libBIE.a

# Dependency
depend .depend : 
	$(CXX) $(CXXFLAGS) -MM $(wildcard *.cc) > .depend

include .depend
