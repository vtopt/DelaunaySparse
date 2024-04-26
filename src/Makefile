FORT = gfortran
CFLAGS = -c
OPTS = -fopenmp 

all: samples samplep test_install
	./test_install

test_install: test_install.f90 delsparse.o bqpd.o lapack.o blas.o
	$(FORT) $(OPTS) test_install.f90 delsparse.o bqpd.o lapack.o blas.o -o test_install

samples: samples.f90 delsparse.o bqpd.o lapack.o blas.o
	$(FORT) $(OPTS) samples.f90 delsparse.o bqpd.o lapack.o blas.o -o samples

samplep: samplep.f90 delsparse.o bqpd.o lapack.o blas.o
	$(FORT) $(OPTS) samplep.f90 delsparse.o bqpd.o lapack.o blas.o -o samplep

delsparse.o: delsparse.f90 
	$(FORT) $(CFLAGS) $(OPTS) delsparse.f90 -o delsparse.o

bqpd.o : bqpd_min/bqpd.f
	$(FORT) $(CFLAGS) $(OPTS) bqpd_min/bqpd.f -o bqpd.o

lapack.o : lapack.f
	$(FORT) $(CFLAGS) $(OPTS) lapack.f -o lapack.o

blas.o : blas.f
	$(FORT) $(CFLAGS) $(OPTS) blas.f -o blas.o

clean:
	rm -f *.o *.mod samples samplep test_install
