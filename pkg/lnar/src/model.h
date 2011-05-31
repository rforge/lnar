#include <gsl/gsl_vector.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#define printf Rprintf

#ifndef STRSPAR
#define STRSPAR
/* Implementation with pointer manipulation */
typedef  void (* funderivs) (int neq, double * t, double *y, double * fout, double *thetas);

typedef struct spar_ {
  double * dataset;
  int   nrows; /* Number of data points*/
  int   ncols; /* Number of Species + 1 (Time) */
  int nmols;   /* System size, e.g. #molecules */
  double * k ; /* Various constants used in the system  */
  int nk; /* Number of constants */
  int odesize; /* Number of elements in the vector ODE */
  int sysdim;  /*The system dimension*/
  int maxiter; /*Hard limit of iterations*/
  double simplexsize; /*Specifies the convergence criterion*/
  double hessianh;/*interval for evaluating the Hessian 
			  matrix (finite differences) */
  double relerr; /*integrator step size: relative tolerance*/
  double abserr; /*integrator step size: absolute tolerance- */

  gsl_vector * x; /* initial point*/
  gsl_vector *ss;  /*  simplex step size */
  double * thetas; /*Current thetas, used in functions*/
  funderivs f;
  int (*jac) (int neq, double t, const double y[], double *dfdy, 
	               double dfdt[], double *params);
  int logc; /* use barrier */
  int method; /* method - concentrations */
  const int *order; /* indicate the order of the functions */
  double logbarp;
  double *thetatrue;
  /*
   * Profile Likelihood
   */
  int pindex; /*index for fixed parameter */
  int pindex2; /*index for fixed 2nd parameter */
  double proftheta; /* Value for the parameter  */
  double proftheta2; /* Value for the 2nd parameter  */
  int verbose;
  int printci; /*Print Confidence Intervals*/
} spar;

#endif

//int derivs(double t, const double * y, double * f,void *params);
//void derivs2(double t, const double *y, double *f, void *params);
int init_model(spar ** mypar, char * filename);
int fin_model( spar ** mypar);
