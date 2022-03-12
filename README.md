# ACM TOMS Algorithm 1012: DELAUNAYSPARSE -- Interpolation via a Sparse Subset of the Delaunay Triangulation

The package DELAUNAYSPARSE (ACM TOMS Algorithm 1012) contains serial and
parallel codes, written in FORTRAN 2003 with OpenMP 4.5, for performing
interpolation in medium to high dimensions via a sparse subset of the
Delaunay triangulation. In addition to the original Fortran source code,
this repository contains a wrapper for Python 3.6+ and C bindings.
Command line drivers are also provided with the original Fortran code.

## Usage

`DELAUNAYSPARSE` contains several modes of operation.

At the most basic level, the two driver subroutines are as follows.
 * `DELAUNAYSPARSES` runs the serial driver to identify the vertices
   of the simplex/simplices containing one or more interpolation points.
   Can also (optionally) be set to compute and return the value of the
   Delaunay interpolant.
 * `DELAUNAYSPARSEP` runs the parallel driver to identify the vertices
   of the simplex/simplices containing one or more interpolation points.
   Can also (optionally) be set to compute and return the value of the
   Delaunay interpolant (must set the `OMP_NUM_THREADS` environment
   variable).

Additionally, two command-line drivers are provided, which read input
from files:
 * `delsparses` (uses the serial driver), and
 * `delsparsep` (uses the parallel driver).

In the `extras` directory, there are two additional interfaces for calling
from C/C++ (`extras/c_binding`) and Python 3 (`extras/delsparsepy`).

Further detailed user information is documented in the USAGE document.

## Organization

The physical organization is as follows.

 * `toms1012` contains the original unmodified Fortran source code, as
   published in ACM TOMS Algorithm 1012. This includes 2 command line drivers
   `samples.f90` (serial driver) and `samplep.f90` (parallel driver), which
   can be used on formatted data files from the command line.
   Comments at the top of each subroutine document their usage.
   See this directory's internal README for further information on
   building, testing, and usage. The directory can be independently
   downloaded and is fully portable.
 * `src` contains the latest project source, which has been configured
   to easily install using modern package managers.
 * `test` contains basic regression test cases for each major mode of
   operation, so that the installation can be tested.
 * `data` contains several real-world data files (exhibiting degeneracy),
   for testing the installation and providing a sample for the CL interface.
 * `extras` contains the C bindings (`c_binding`) and DelaunaySparse for
   Python (`delsparsepy`). For convenience, copies of all source code and
   dependencies are duplicated in each of these directories.
 * A GNU Makefile is provided for building, installing, and running tests.
   Define your Fortran compiler and options at the top of the file, and
   (optionally) set the install directory by setting the $(PREFIX) variable.
   Then use `make` to build binaries, `make install` to install binaries,
   and `make test_install` to test the installation.
   Binary executables will install in `$(PREFIX)/bin`, headers will install
   to `$(PREFIX)/include`, and shared libraries will install to
   `$(PREFIX)/lib`. Running these commands also builds/installs/tests the
   C interface, `delsparsec`. The Python extras must be installed separately.
 * USAGE provides additional detailed user information.

## Citation

To cite this work, please use

```
@article{TOMSalg1012,
author = {Chang, Tyler H. and Watson, Layne T. and Lux, Thomas C. H. and Butt, Ali R. and Cameron, Kirk W. and Hong, Yili},
title = {Algorithm 1012: {DELAUNAYSPARSE}: {I}nterpolation via a Sparse Subset of the {D}elaunay Triangulation in Medium to High Dimensions},
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
```

## Inquiries

For further inquiries, contact
Tyler Chang, tchang@anl.gov.
