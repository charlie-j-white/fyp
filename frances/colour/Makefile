# compiles the c programs so they can be called from fortran


CC = gcc
CFLAGS = -lm -O3 -Wall

FC = gfortran
FFLAGS = -O3



%.o: %.c
	$(CC) -c $< $(CFLAGS)


%.o: %.f
	$(FC) -ffree-form -c $< $(FFLAGS)

OBJS = main.o\
       colour.o\
       lodepng.o




main: $(OBJS)
	$(FC) -o $@ $^ $(FFLAGS)









