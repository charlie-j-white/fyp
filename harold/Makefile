#
SHELL=/bin/sh
#
F77 = gfortran
#
#F77FLAGS = -O3 -C -Minform=inform
F77FLAGS = -O3  
#
.SUFFIXES:	.o .f

.f.o:;	$(F77) -c $(F77FLAGS) $<

SRCS=	boundaries.f\
	griddata.f\
	initialise.f\
	nozzle_main.f\
	readconfig.f\
	resid.f\
	stagedata.f\
	timestep.f\
	update.f\
	wrapper.f\
	cost.f\
	costMERGE_d-all.f\
	jRES.f\
	resMERGE_d.f

OBJS=	boundaries.o\
	griddata.o\
	initialise.o\
	nozzle_main.o\
	readconfig.o\
	resid.o\
	stagedata.o\
	timestep.o\
	update.o\
	wrapper.o\
	cost.o\
	costMERGE_d-all.o\
	jRES.o\
	resMERGE_d.o

EXEC =  nozzle.x

all:	$(OBJS)
	$(F77) $(F77FLAGS) -o $(EXEC) $(OBJS) $(LFLAGS)

