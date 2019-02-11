
FC = gfortran
FFLAGS = -O3

.f.o:
	$(FC) $(FFLAGS) -c $<

OBJS= 	main.o\
	meshing.o\
	initialise.o\
	boundaries.o\
	update.o\
	fvector.o\
	qvector.o

main: $(OBJS)
	$(FC) -o $@ $^ $(FFLAGS)

.PHONY: clean

clean: 
	rm $(OBJS)



