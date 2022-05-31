#ifndef DELSPARSEC
#define DELSPARSEC

// serial subroutine: no optional arguments
extern void c_delaunaysparses(int *d, int *n, double pts[], int *m, double q[],
                              int simps[], double weights[], int ierr[]);

// serial: compute interpolant values
extern void c_delaunaysparses_interp(int *d, int *n, double pts[], int *m,
                              double q[], int simps[], double weights[],
                              int ierr[], int *ir, double interp_in[],
                              double interp_out[]);

// serial: optional arguments, no interpolant values
extern void c_delaunaysparses_opts(int *d, int *n, double pts[], int *m,
                                   double q[],int simps[], double weights[],
                                   int ierr[], double *eps, double *extrap,
                                   double rnorm[], int *ibudget, bool *chain,
                                   bool *exact);

// serial: optional arguments and compute interpolant values
extern void c_delaunaysparses_interp_opts(int *d, int *n, double pts[], int *m,
                                          double q[],int simps[],
                                          double weights[], int ierr[],
                                          int *ir, double interp_in[],
                                          double interp_out[], double *eps,
                                          double *extrap, double rnorm[],
                                          int *ibudget, bool *chain,
                                          bool *exact);


// parallel: no optional arguments
extern void c_delaunaysparsep(int *d, int *n, double pts[], int *m, double q[],
                              int simps[], double weights[], int ierr[]);

// parallel: compute interpolant values
extern void c_delaunaysparsep_interp(int *d, int *n, double pts[], int *m,
                              double q[], int simps[], double weights[],
                              int ierr[], int *ir, double interp_in[],
                              double interp_out[]);

// parallel: optional arguments, no interpolant values
extern void c_delaunaysparsep_opts(int *d, int *n, double pts[], int *m,
                                   double q[],int simps[], double weights[],
                                   int ierr[], double *eps, double *extrap,
                                   double rnorm[], int *ibudget, bool *chain,
                                   bool *exact, int *pmode);

// parallel: optional arguments and compute interpolant values
extern void c_delaunaysparsep_interp_opts(int *d, int *n, double pts[], int *m,
                                          double q[],int simps[],
                                          double weights[], int ierr[],
                                          int *ir, double interp_in[],
                                          double interp_out[], double *eps,
                                          double *extrap, double rnorm[],
                                          int *ibudget, bool *chain,
                                          bool *exact, int *pmode);

#endif
