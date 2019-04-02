#**********************************************************************#
#                          Charlie Anderson                            #
#                  University of Bristol, March 2019                   #
#**********************************************************************#


# define compiler options
FC =gfortran
FLIBS =
FFLAGS =-O3 -Wall


# name which object files to use
_OBJS = wrapper.o\
       main.o\
       meshing.o\
       split_fwd.o\
       split_rev.o\
       debug_cell.o\
       postprocess.o\
       initialise.o\
       resid.o\
       boundaries.o\
       update.o\
       pressure.o\
       fg_vector.o\
       itinfo.o\
       jst_calcs.o\
       timestep.o


# specify directory for object files and put them there
ODIR =OBJ
OBJS = $(patsubst %,$(ODIR)/%,$(_OBJS))


# first rule, compile all *.f files
$(ODIR)/%.o: %.f
	$(FC) -c -o $@ $< $(FFLAGS)


# main rule, performs linking
main: $(OBJS)
	$(FC) -o $@ $^ $(FFLAGS) $(FLIBS)


# last rule, `refreshes' the build by deleting the necessary files
.PHONY: clean
clean:
	rm $(OBJS) main
