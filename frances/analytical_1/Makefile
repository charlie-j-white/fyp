
FC = gfortran
FFLAGS = -O3

.f.o:
	$(FC) $(FFLAGS) -c $<

OBJS= 	main.o\
	update.o

main: $(OBJS)
	$(FC) -o $@ $^ $(FFLAGS)

.PHONY: clean

clean: 
	rm $(OBJS)



