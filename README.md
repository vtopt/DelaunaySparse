# ACM TOMS Algorithm 1012: DELAUNAYSPARSE -- Interpolation via a Sparse Subset of the Delaunay Triangulation

The package DELAUNAYSPARSE (ACM TOMS Algorithm 1012) contains serial and
parallel codes, written in FORTRAN 2003 with OpenMP 4.5, for performing
interpolation in medium to high dimensions via a sparse subset of the
Delaunay triangulation. In addition to the original Fortran source code,
this repository contains a wrapper for Python 3.6+ and C bindings.
Command line drivers are also provided with the original Fortran code.

## Organization and Usage

The physical organization is as follows. Note that each of the following
directories could be independently downloaded.

 * `src` contains the original unmodified Fortran source code, as published
   in ACM TOMS Algorithm 1012. This includes 2 command line drivers
   `samples.f90` (serial driver) and `samplep.f90` (parallel driver), which
   can be used on formatted data files from the command line.
   Comments at the top of each subroutine document their usage.
   See this directory's internal README for further information on
   building, testing, and usage.
 * `python` contains a Python3 wrapper for the Fortran code, allowing
   DELAUNAYSPARSE to be directly imported as a Python package. This wrapper
   was created by modifying the output generated by fmodpy. The script
   `example.py` demonstrates its usage. For convenience, copies of all
   Fortran code that is used by the Python wrapper are also included in
   this directory.
 * `c_binding` contains C bindings for several variations of the main
   Fortran subroutines, as well as copies of the Fortran source code.
   A test file `test_install.c` can be used for usage examples. This
   directory's internal README also contains best practices when calling
   Fortran from C/C++.

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
