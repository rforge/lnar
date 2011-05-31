/***************************************************************************************
 *  Multivariate Normal log-density function
 *  Using GSL -> www.gnu.org/software/gsl
 *
 *  Copyright (C) 2006  Ralph dos Santos Silva
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 *  AUTHOR 
 *     Ralph dos Santos Silva,  [EMAIL PROTECTED]
 *     March, 2006       
 *     Vasileios Giagos: a minor modification to return the log-density
***************************************************************************************/

//#define PRINTSTIFF
#include "mnormal.h"
#include <gsl/gsl_math.h>
#include "rand.h"
#ifdef PRINTSTIFF
#include <gsl/gsl_eigen.h>
#endif

double ldmvnorm(const int n, const gsl_vector *x, const gsl_vector *mean, const gsl_matrix *var){
    /* log-multivariate normal density function    */
    /*
     *	n	dimension of the random vetor
     *	mean	vector of means of size n
     *	var	variance matrix of dimension n x n
     */
    int s;
    double ax,ay;
    gsl_vector *ym, *xm;
    gsl_matrix *work = gsl_matrix_alloc(n,n), 
	                  *winv = gsl_matrix_alloc(n,n);
#ifdef PRINTSTIFF
    /* Print Stiffness indicator S=max(eigen)/min(eigen)*/
    gsl_vector *eval = gsl_vector_alloc (n);
    gsl_matrix *evec = gsl_matrix_alloc (n, n);
    gsl_matrix_memcpy(work,var);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(n);
    gsl_eigen_symmv (work, eval, evec, w);
    gsl_eigen_symmv_free (w);
    gsl_eigen_symmv_sort (eval, evec, 
	    GSL_EIGEN_SORT_ABS_ASC);
    printf("%f ",
	    gsl_vector_get(eval,n-1) / gsl_vector_get(eval,0));
    gsl_vector_free(eval);
    gsl_matrix_free(evec);
#endif

    gsl_permutation *p = gsl_permutation_alloc(n);
    gsl_matrix_memcpy( work, var );
    gsl_linalg_LU_decomp( work, p, &s );
    gsl_linalg_LU_invert( work, p, winv );
    ax = gsl_linalg_LU_det( work, s );
    gsl_matrix_free( work );
    gsl_permutation_free( p );

    xm = gsl_vector_alloc(n);
    gsl_vector_memcpy( xm, x);
    gsl_vector_sub( xm, mean );
    ym = gsl_vector_alloc(n);
    gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
    gsl_matrix_free( winv );
    gsl_blas_ddot( xm, ym, &ay);
    gsl_vector_free(xm);
    gsl_vector_free(ym);
    /* 
     * ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),n)*ax );
     */
    ay = -0.5*( ay + n*log(2*M_PI) + log(ax) );
    return ay;
}

int rmvnorm(const gsl_rng *r, const int n, const gsl_vector *mean, 
		const gsl_matrix *var, gsl_vector *result){
    /* multivariate normal distribution random number generator */
    /*
     *	n	dimension of the random vetor
     *	mean	vector of means of size n
     *	var	variance matrix of dimension n x n
     *	result	output variable with a sigle random vector normal distribution generation
     */
    int k;
    gsl_matrix *work = gsl_matrix_alloc(n,n);



    gsl_matrix_memcpy(work,var);
    gsl_linalg_cholesky_decomp(work);

    for(k=0; k<n; k++)
	gsl_vector_set( result, k, gsl_ran_ugaussian(r) );

    gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, work, result );
    gsl_vector_add(result,mean);
    gsl_matrix_free(work);
    return 0;
}
