
Simulation += Simulation_data.o profile_helmholtz.o nr.o nrutil.o nrtype.o newt.o fdjac.o ludcmp.o lubksb.o fmin.o lnsrch.o sim_newt_functions.o

Gravity += Gravity_sendOutputData.o

Simulation_init.o : Simulation_data.o sim_newt_functions.o
Simulation_initBlock.o : Simulation_data.o

nrutil.o : nrtype.o
newt.o : nr.o nrutil.o fdjac.o
fmin.o : nrutil.o
fdjac.o : nrutil.o
ludcmp.o : nrutil.o
lnsrch.o : nrutil.o
