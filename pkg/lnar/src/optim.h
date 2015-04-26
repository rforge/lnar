#include "lik.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

double my_f (const gsl_vector *v, void *params);
double my_f2 (const gsl_vector *v, void *params);
SEXP varcov(spar **mypar, const gsl_vector *v ,int printh);
SEXP genres(spar **mypar, const gsl_vector *v ,int printh);
int optim(spar **mypar);
int optim2(spar **mypar);
double proflik(spar ** mypar);
double proflik2(spar ** mypar);
