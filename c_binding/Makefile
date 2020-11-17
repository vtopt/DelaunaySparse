FORT = gfortran
CC = gcc
CFLAGS = -c
OPTS = -fopenmp 
LEGACY = -std=legacy

all: test_install.o delsparse_bind_c.o delsparse.o slatec.o lapack.o blas.o
	$(FORT) $(OPTS) test_install.o delsparse_bind_c.o delsparse.o slatec.o lapack.o blas.o -o test_install
	./test_install

test_install.o: test_install.c
	$(CC) $(CFLAGS) $(OPTS) test_install.c -o test_install.o

delsparse_bind_c.o: delsparse_bind_c.f90 delsparse.o
	$(FORT) $(CFLAGS) $(OPTS) delsparse_bind_c.f90 -o delsparse_bind_c.o

delsparse.o: delsparse.f90 
	$(FORT) $(CFLAGS) $(OPTS) delsparse.f90 -o delsparse.o

slatec.o : slatec.f
	$(FORT) $(CFLAGS) $(OPTS) $(LEGACY) slatec.f -o slatec.o

lapack.o : lapack.f
	$(FORT) $(CFLAGS) $(OPTS) lapack.f -o lapack.o

blas.o : blas.f
	$(FORT) $(CFLAGS) $(OPTS) blas.f -o blas.o

clean:
	rm -f *.o *.mod test_install
