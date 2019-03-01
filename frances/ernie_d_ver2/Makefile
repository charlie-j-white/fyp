
FC = gfortran
#FFLAGS = -O3 --std=legacy -Wall
FFLAGS = -O3

.f.o:
	$(FC) $(FFLAGS) -c $<

#OBJS:= $(wildcard *.o)
OBJS= 	main.o\
	meshing.o\
	initialise.o\
	boundaries.o\
	update.o\
	fvector.o\
	qvector.o\
	jcost.o\
	wrapper.o\
	MERGED_d-all.o

main: $(OBJS)
	$(FC) -o $@ $^ $(FFLAGS)

.PHONY: clean

clean: 
	rm main solution.plt MERGED.f $(OBJS)



