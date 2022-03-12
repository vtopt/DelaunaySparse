# Your Fortran 2003 compliant compiler. Must support C bindings to build
# delsparsec library (in extras/c_binding).
FORT = gfortran

# Set the prefix for your install directory. Installs in-place by default.
PREFIX = $(PWD)

# Compiler flags defined below.

# Build without linking
CFLAGS = -c

# $OPTS must contain the flag for building OpenMP threadsafe (-fopenmp on GNU).
# You can add additional flags (such as -O3) for code optimization
OPTS = -fopenmp -O3 -fPIC

# Link shared objects
LFLAGS = -shared

# Legacy flag, used for suppressing warnings when building SLATEC
LEGACY = -std=legacy

# Dependencies: Replace with appropriate linker flag (e.g., -llapack -lblas)
# if you have these libraries already installed on your computer.
# Otherwise, the value below will use the included minimal copies, taken
# from the public domain.
# Note that there is a known issue that occurs during extrapolation when
# linking against the public version of lapack.f.
LIBS = dependencies/lapack.f dependencies/blas.f


# Build formula:

# Make all library/module files

all: src/libdelsparse.so src/delsparses src/delsparsep extras/c_binding/libdelsparsec.so

# Copy all libraries, module files, and binaries into their install locations

install: include bin lib
	cp src/libdelsparse.so $(PREFIX)/lib/libdelsparse.so
	cp src/delsparse_mod.mod $(PREFIX)/include/delsparse_mod.mod
	cp src/real_precision.mod $(PREFIX)/include/real_precision.mod
	cp src/delsparses $(PREFIX)/bin/delsparses
	cp src/delsparsep $(PREFIX)/bin/delsparsep
	cp extras/c_binding/libdelsparsec.so $(PREFIX)/lib/libdelsparsec.so
	cp extras/c_binding/delsparse.h $(PREFIX)/include/delsparse.h

# Test installation

test_install: test/test_install.f90 test/test_c_install.c include/delsparse_mod.mod include/delsparse.h lib/libdelsparse.so lib/libdelsparsec.so bin/delsparses bin/delsparsep
	$(FORT) $(OPTS) test/test_install.f90 -I$(PREFIX)/include -L$(PREFIX)/lib -ldelsparse -o test/test_install
	export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(PREFIX)/lib && test/test_install
	$(FORT) $(OPTS) test/test_c_install.c -I$(PREFIX)/include -L$(PREFIX)/lib -ldelsparsec -o test/test_c_install
	export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(PREFIX)/lib && test/test_c_install
	test/test_bin.sh

# Make shared libs

extras/c_binding/libdelsparsec.so: src/dependencies/slatec.o extras/c_binding/delsparsec.o
	cp src/dependencies/slatec.o extras/c_binding/dependencies/slatec.o
	cd extras/c_binding && $(FORT) $(LFLAGS) $(OPTS) delsparsec.o delsparse.o dependencies/slatec.o $(LIBS) -o libdelsparsec.so

src/libdelsparse.so: src/dependencies/slatec.o src/delsparse.o
	cd src && $(FORT) $(LFLAGS) $(OPTS) delsparse.o dependencies/slatec.o $(LIBS) -o libdelsparse.so

# Make bin execs

src/delsparses: src/samples.f90 src/delsparse.o src/dependencies/slatec.o
	cd src && $(FORT) $(OPTS) samples.f90 delsparse.o dependencies/slatec.o $(LIBS) -o delsparses

src/delsparsep: src/samplep.f90 src/delsparse.o src/dependencies/slatec.o
	cd src && $(FORT) $(OPTS) samplep.f90 delsparse.o dependencies/slatec.o $(LIBS) -o delsparsep

# Make C bindings

extras/c_binding/delsparsec.o: extras/c_binding/delsparse_bind_c.f90
	cd extras/c_binding && $(FORT) $(CFLAGS) $(OPTS) delsparse.f90 -o delsparse.o
	cd extras/c_binding && $(FORT) $(CFLAGS) $(OPTS) delsparse_bind_c.f90 -o delsparsec.o

# Make delsparse.o and slatec.o

src/delsparse.o: src/delsparse.f90 
	cd src && $(FORT) $(CFLAGS) $(OPTS) delsparse.f90 -o delsparse.o

src/dependencies/slatec.o : src/dependencies/slatec.f
	cd src/dependencies && $(FORT) $(CFLAGS) $(OPTS) $(LEGACY) slatec.f -o slatec.o

# Make install directories

include:
	mkdir include

lib:
	mkdir lib

bin:
	mkdir bin

# Clean command

clean:
	cd src && rm -f *.o *.mod *.so
	cd src/dependencies && rm -f *.o
	cd extras/c_binding && rm -f *.o *.mod *.so
	cd extras/c_binding/dependencies && rm -f *.o
	cd lib && rm -f *.so
	cd bin && rm -f delsparses delsparsep
	cd include && rm -f delsparse.h *.mod
	cd test && rm -f test_install test_c_install
