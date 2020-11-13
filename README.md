# DELAUNAYSPARSE: Interpolation via the Delaunay Triangulation.

The package DELAUNAYSPARSE contains serial and parallel codes, written
in FORTRAN 2003 with OpenMP, for performing interpolation in medium to
high dimensions via a sparse subset of the Delaunay triangulation. The
serial driver subroutine is `DELAUNAYSPARSES` and the parallel driver is
`DELAUNAYSPARSEP`. Both subroutines use the `REAL_PRECISION` module from
HOMPACK90 (ACM TOMS Algorithm 777) for approximately 64-bit precision
on all known machines, and the SLATEC subroutine `DWNNLS` (ACM TOMS
Algorithm 587) for solving an inequality constrained least squares
problem. Additionally, `DELAUNAYSPARSE` depends on several BLAS and LAPACK
subroutines. The module `DELSPARSE_MOD` contains the `REAL_PRECISION` (R8)
data type, and interface blocks for `DELAUNAYSPARSES`, `DELAUNAYSPARSEP`,
and `DWNNLS`. Comments at the top of each subroutine document their
usage, and examples demonstrating their usage are provided by the
sample command line programs in `samples.f90` and `samplep.f90`.

The physical organization is as follows:

 * The file `deslparse.f90` contains the module `REAL_PRECISION`,
   `DELSPARSE_MOD`, and the driver subroutines `DELAUNAYSPARSES`, and
   `DELAUNAYSPARSEP`.
 * The file `slatec.f` contains the subroutine `DWNNLS` and its dependencies
   from the SLATEC library. This library has been slightly modified to
   comply with the modern Fortran standards. Additionally, legacy
   implementations of the BLAS subroutines `DROTM` and `DTROMG` have been
   included under different names to avoid dependency issues.
 * The file `samples.f90` contains a sample main program demonstrating the
   usage of `DELAUNAYSPARSES`, with optional arguments.
 * The file `samplep.f90` contains a sample main program demonstrating the
   usage of `DELAUNAYSPARSEP`, with optional arguments.
 * The file `test_install.f90` contains a simple test program that checks
   whether the installation of `DELAUNAYSPARSE` appears correct, based
   on the output to a small interpolation/extrapolation problem.
 * The file `sample_input2d.dat` contains a sample 2-dimensional input
   data set for `samples.f90` and `samplep.f90`.
 * The file `sample_input4d.dat` contains a sample 4-dimensional input
   data set for `samples.f90` and `samplep.f90`.
 * The files `lapack.f` and `blas.f` contain all LAPACK and BLAS
   subroutines that are referenced (both directly and indirectly) in
   `DELAUNAYSPARSE`.
 * A sample GNU `Makefile` is provided.

From here on, the files `samples.f90` and `samplep.f90` will be referred
to collectively as `sample{s|p}.f90` and the files `sample_input2d.dat`
and `sample_input4d.dat` will be referred to collectively as
`sample_input{2|4}d.dat`.

To check that the installation of `DELAUNAYSPARSES` and `DELAUNAYSPARSEP` is
correct, assuming that your Fortran compiler allows mixing fixed format
`.f` and free format `.f90` files in the same compile command, use the command

``
$FORT $OPTS delsparse.f90 slatec.f lapack.f blas.f test_install.f90 \
  -o test_install $LIBS
``

where `$FORT` is a Fortran 2003 compliant compiler supporting OpenMP
4.5, `$OPTS` is a list of compiler options, and `$LIBS` is a list of
flags to link the BLAS and LAPACK libraries, if those exist on your
system (in which case the files `blas.f` and `lapack.f` can be omitted
from the compile command). To run the parallel code, `$OPTS` must
include the compiler option for OpenMP.

Then run the tests using

``
./test_install
``

To compile and link the sample main programs `sample{s|p}.f90`, use

``
$FORT $OPTS delsparse.f90 slatec.f lapack.f blas.f sample{s|p}.f90 \
  -o sample{s|p} $LIBS
``

similar to above.  To run a sample main program, use

``
./sample{s|p} sample_input{2|4}d.dat
``

where `sample_input{2|4}d.dat` could be replaced by any other similarly
formatted data file.

---------------------------------------------------------------------------

For further inquiries, contact
Tyler Chang, tchang@anl.gov.

---------------------------------------------------------------------------

To cite this work, please use

``
@article{TOMSalg1012,
author = {Chang, Tyler H. and Watson, Layne T. and Lux, Thomas C. H.
and Butt, Ali R. and Cameron, Kirk W. and Hong, Yili},
title = {Algorithm 1012: DELAUNAYSPARSE: Interpolation via a Sparse Subset of
the Delaunay Triangulation in Medium to High Dimensions},
year = {2020},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {46},
number = {4},
doi = {10.1145/3422818},
journal = {ACM Trans. Math. Softw.},
month = {nov},
articleno = {38},
numpages = {20}
}
``
