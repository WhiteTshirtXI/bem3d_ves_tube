This code is a modified version of Hong Zhao's blood code for large-scale
simulation of the microcirculation. It is reformatted to model the
low-Reynolds number motion of a vesicle through a channel.

In short, all references to the "Cell" class are changed to "Vesicle" class.
The stokesys_flowsolver.cc code is still used to assemble the resistance
equation, but references to "Rigid" are not utilized. A single vesicle is
placed in the channel.

Changes made (7/25/2014 by Andrew Spann):
vessim/channel.cc
lib/stokesys.cc
include/stokesys.h
lib/stokesys_flowsolver.cc
vessim/vesinitcond.cc

Things to change: stokesys_flowsolver.cc
- this library assembles the resistance matrix corresponding to cell-cell,
  cell-rigid, cell-wall, etc. interactions
- look at solveFlow(), which will call calclhs and calcrhs (LHS and RHS of the
  resistance equation)
- need to get rid of anything that refers to an elastic force
- need to add in things that solve for surface tension using vesicle functions
- don't need to change vesicle.cc or vesicle.h, but will be calling things in
  those functions
- vesys.cc is the template you'll look at in order to make the appropriate
  changes to stokesys_flowsolver.cc (try to replicate the content of vesys.cc
into stokesys_flowsolver)
- goal is to make the content of the flow solver in vesys.cc (e.g. the
  solveVel function) play nice with the walls
