FORT = gfortran
CC = gcc
CFLAGS = -c
OPTS = -fopenmp 
LIBS = dependencies/blas.f dependencies/lapack.f
LEGACY = -std=legacy

all: test_install.o delsparse_bind_c.o delsparse.o dependencies/slatec.o delsparse.h
	$(FORT) $(OPTS) test_install.o delsparse_bind_c.o delsparse.o dependencies/slatec.o $(LIBS) -o test_install
	./test_install

test_install.o: test_install.c
	$(CC) $(CFLAGS) $(OPTS) test_install.c -o test_install.o

delsparse_bind_c.o: delsparse_bind_c.f90 delsparse.o
	$(FORT) $(CFLAGS) $(OPTS) delsparse_bind_c.f90 -o delsparse_bind_c.o

delsparse.o: delsparse.f90 
	$(FORT) $(CFLAGS) $(OPTS) delsparse.f90 -o delsparse.o

dependencies/slatec.o : dependencies/slatec.f
	cd dependencies && $(FORT) $(CFLAGS) $(OPTS) $(LEGACY) slatec.f -o slatec.o

clean:
	rm -f *.o *.mod test_install
