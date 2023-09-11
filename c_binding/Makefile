FORT = gfortran
CC = gcc
CFLAGS = -c
OPTS = -fopenmp 
LIBS = blas.f lapack.f

all: test_install.o delsparse_bind_c.o delsparse.o bqpd.o delsparse.h
	$(FORT) $(OPTS) test_install.o delsparse_bind_c.o delsparse.o bqpd.o $(LIBS) -o test_install
	./test_install

test_install.o: test_install.c
	$(CC) $(CFLAGS) test_install.c -o test_install.o

delsparse_bind_c.o: delsparse_bind_c.f90 delsparse.o
	$(FORT) $(CFLAGS) $(OPTS) delsparse_bind_c.f90 -o delsparse_bind_c.o

delsparse.o: delsparse.f90 
	$(FORT) $(CFLAGS) $(OPTS) delsparse.f90 -o delsparse.o

bqpd.o : bqpd_min/bqpd.f
	$(FORT) $(CFLAGS) $(OPTS) $(LEGACY) bqpd_min/bqpd.f -o bqpd.o

clean:
	rm -f *.o *.mod test_install
