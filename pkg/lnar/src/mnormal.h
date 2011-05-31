#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "rand.h"

double ldmvnorm(const int n, const gsl_vector *x, 
	const gsl_vector *mean, const gsl_matrix *var);
int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, 
		const gsl_matrix *var, gsl_vector * res);
