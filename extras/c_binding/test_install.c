#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "delsparse.h"

int main() {
   // Set the problem dimensions
   int n = 50, d = 5, m = 10, ir = 2;

   // Generate random data in the unit cube
   double data[n*d];
   for (int i = 0; i < n*d; i++)
      data[i] = rand();

   // Generate interpolation points
   double interp[m*d];
   for (int i = 0; i < m*d; i++)
      interp[i] = 0.25 + 0.5 * rand();

   // Generate response values
   double interp_in[n*ir];
   for (int i = 0; i < n*ir; i++)
      interp_in[i] = rand();

   // Allocate the output arrays
   int simps[m*(d+1)], ierr[m];
   double weights[m*(d+1)], interp_out[m*ir], rnorm[m];

   // Set the optional input parameters
   bool chain = false, exact = true;
   int ibudget = 10000, pmode = 1;
   double eps = 0.00000001, extrap = 0.1;

   // Call the serial C interface with no options
   c_delaunaysparses(&d, &n, data, &m, interp, simps, weights, ierr);

   // Check for errors
   for (int i = 0; i < m; i++) {
      if (ierr[i] > 2) {
         printf("Error %i occurred while testing c_delaunaysparses"
                " with no optional arguments\n\n",
                ierr[i]);
         return -1;
      }
   }

   // Call the serial C interface and compute interpolant values
   c_delaunaysparses_interp(&d, &n, data, &m, interp, simps, weights, ierr,
                            &ir, interp_in, interp_out);

   // Check for errors
   for (int i = 0; i < m; i++) {
      if (ierr[i] > 2) {
         printf("Error %i occurred while testing c_delaunaysparses"
                " and computing interpolant values\n\n", ierr[i]);
         return -1;
      }
   }

   // Call the serial C interface with optional inputs
   c_delaunaysparses_opts(&d, &n, data, &m, interp, simps, weights, ierr,
                          &eps, &extrap, rnorm, &ibudget, &chain, &exact);

   // Check for errors
   for (int i = 0; i < m; i++) {
      if (ierr[i] > 2) {
         printf("Error %i occurred while testing c_delaunaysparses"
                " with optional arguments\n\n", ierr[i]);
         return -1;
      }
   }

   // Call the serial C interface with optional inputs and interpolation
   c_delaunaysparses_interp_opts(&d, &n, data, &m, interp, simps, weights,
                                 ierr, &ir, interp_in, interp_out, &eps,
                                 &extrap, rnorm, &ibudget, &chain, &exact);

   // Check for errors
   for (int i = 0; i < m; i++) {
      if (ierr[i] > 2) {
         printf("Error %i occurred while testing c_delaunaysparses"
                " with optional arguments and computing the interpolant\n\n",
                ierr[i]);
         return -1;
      }
   }


   // Call the parallel C interface with no options
   c_delaunaysparsep(&d, &n, data, &m, interp, simps, weights, ierr);

   // Check for errors
   for (int i = 0; i < m; i++) {
      if (ierr[i] > 2) {
         printf("Error %i occurred while testing c_delaunaysparsep"
                " with no optional arguments\n\n",
                ierr[i]);
         return -1;
      }
   }

   // Call the parallel C interface and compute interpolant values
   c_delaunaysparsep_interp(&d, &n, data, &m, interp, simps, weights, ierr,
                            &ir, interp_in, interp_out);

   // Check for errors
   for (int i = 0; i < m; i++) {
      if (ierr[i] > 2) {
         printf("Error %i occurred while testing c_delaunaysparsep"
                " and computing interpolant values\n\n", ierr[i]);
         return -1;
      }
   }

   // Call the parallel C interface with optional inputs
   c_delaunaysparsep_opts(&d, &n, data, &m, interp, simps, weights, ierr,
                          &eps, &extrap, rnorm, &ibudget, &chain, &exact,
                          &pmode);

   // Check for errors
   for (int i = 0; i < m; i++) {
      if (ierr[i] > 2) {
         printf("Error %i occurred while testing c_delaunaysparsep"
                " with optional arguments\n\n", ierr[i]);
         return -1;
      }
   }

   // Call the parallel C interface with optional inputs and interpolation
   c_delaunaysparsep_interp_opts(&d, &n, data, &m, interp, simps, weights,
                                 ierr, &ir, interp_in, interp_out, &eps,
                                 &extrap, rnorm, &ibudget, &chain, &exact,
                                 &pmode);

   // Check for errors
   for (int i = 0; i < m; i++) {
      if (ierr[i] > 2) {
         printf("Error %i occurred while testing c_delaunaysparsep"
                " with optional arguments and computing the interpolant\n\n",
                ierr[i]);
         return -1;
      }
   }


   // If we made it this far, the build was successful
   printf("The build appears to be successful\n\n");
   return 0;
}
