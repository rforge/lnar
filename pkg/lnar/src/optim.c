#include "optim.h"
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "hessian.h"

//#define SIMPLEX2
//#define DEBUG 1
//#define NORMSCALE2
/*
include <mcheck.h>
*/

#ifdef HAVEOOL
void iteration_echo( ool_conmin_minimizer *M );
#endif
//For the log-scale
double my_f (const gsl_vector *v, void *params){

    spar * myparams= (spar *) params;
    double thetas[myparams->sysdim];
    double scale[myparams->sysdim];
    static double sum;
    int i;
    double lmols=log(myparams->nmols);
    
    sum=0;
    for(i=0;i< myparams->sysdim ;i++)
    {
	thetas[i]=exp(gsl_vector_get(v,i));
	
	if(myparams->logc==1){
	    if(myparams->method==1)
	    {
		//scale[i]=thetas[i]* pow(myparams->nmols,1-myparams->order[i]);
		scale[i]=exp(gsl_vector_get(v,i) -(myparams->order[i]-1)*lmols);
	    }
	    else scale[i]=thetas[i];
          sum+= log(scale[i])+ log(1-scale[i]/3) ;
           //sum+= log(3-scale[i]) ;
	   //sum+= log(3-thetas[i]) ;
        }
	//printf("%4.20f,",thetas[i]);
    }
    //printf("\n");
    return -llik(&myparams,thetas) -myparams->logbarp*sum;
    //return -llik(&myparams,thetas);
}

/*Normal Scale*/
double my_f2 (const gsl_vector *v, void *params){

    spar * myparams= (spar *) params;
    return -llik(&myparams,v->data);
}



double mlik(double x, void *params)
{
  spar* myparams=(spar *) params;
  double thetas[myparams->sysdim];
  int i;
  for(i=0;i<myparams->sysdim;i++)
    {
      thetas[i]=exp(myparams->thetas[i]);
    }
  thetas[myparams->pindex]=exp(x);
  return -llik(&myparams,thetas);  
}

double mlik2(double x, void *params)
{
  spar* myparams=(spar *) params;
  double thetas[myparams->sysdim];
  int i;
  for(i=0;i<myparams->sysdim;i++)
    {
      thetas[i]=exp(myparams->thetas[i]);
    }
  thetas[myparams->pindex]=exp(x);
  return -llik(&myparams,thetas);  
}
/* Numerical Gradient */
 void numgrad (const gsl_vector *x, void * params, gsl_vector *g)
{
  spar* myparams=(spar *) params;
  int i;
  gsl_function F;
  F.function= &mlik;
  F.params= params;
  myparams->thetas=x->data;
  double abserr,result;
  for(i=0;i<myparams->sysdim;i++)
    {
      myparams->pindex=i;
      //gsl_deriv_central(&F,gsl_vector_get(x,i),1e-9,&result,&abserr);
      /*
	for dataset gill-lar-7 (tol =.1) 
	  1e-4: 430 iters
	  1e-5: 399 iters
	  ie-6: 227 iters
      */
      gsl_deriv_central(&F,gsl_vector_get(x,i),1e-5,&result,&abserr);
      //printf("%1.1e ",abserr);
      gsl_vector_set(g,i,result);
    }
  //printf("\n");
}

/* Compute both f and df together. */
 void my_fdf (const gsl_vector *x, void *params, 
	double *f, gsl_vector *df) 
{
  *f = my_f(x, params); 
  numgrad(x, params, df);
}




//For profile likelihood
static double my_pl(const gsl_vector *v, void *params)
{

    spar * myparams= (spar *) params;
    double thetas[myparams->sysdim];
    double scale[myparams->sysdim];
    double sum=0;
    int i,j;
    for(i=0;i< v->size;i++)
    {
	if(myparams->pindex!=i) 
	{
	    thetas[i]=exp(gsl_vector_get(v,i));
        } else 
	{
	    thetas[i]=myparams->proftheta;

	    for(j=(i+1);j<myparams->sysdim;j++)
	    {    
		thetas[j]=exp(gsl_vector_get(v,j-1));
	    }
	    i=v->size;
	}
    }
    if(v->size==myparams->pindex) thetas[myparams->pindex]=myparams->proftheta;
#ifdef DEBUG
    for(i=0;i<myparams->sysdim;i++) printf("%f ",thetas[i]);
    printf("\n");
#endif

    return -llik(&myparams,thetas);
}

//For 2d profile likelihood
static double my_pl2(const gsl_vector *v, void *params)
{

    spar * myparams= (spar *) params;
    double thetas[myparams->sysdim];
    double scale[myparams->sysdim];
    double sum=0;
    register int i,j=0;
    
    for(i=0;i< v->size;i++)
    {
	if(j==myparams->pindex) 
	{
	    thetas[j]=myparams->proftheta;
	    j++;
	}
	if(j==myparams->pindex2)
	{
	    thetas[j]=myparams->proftheta2;
	    j++;
	}
	thetas[j]=exp(gsl_vector_get(v,i));
	j++;
    }
    return -llik(&myparams,thetas);
}

SEXP varcov (spar **mypar, const gsl_vector *v , int printh)
{
    spar * myparams=*mypar;
    int i,j,cover=0;
    double tmp;
    SEXP list,list_names,uconf,est,lconf;
    double * p_uconf,*p_est,*p_lconf;
    char *names[3]={"UC","ES","LC"};

    /*START OF SE EVALUATION*/
    gsl_matrix *m=gsl_matrix_alloc(myparams->sysdim,myparams->sysdim);
    gsl_matrix *minv=gsl_matrix_alloc(myparams->sysdim,myparams->sysdim);


    hessian(v,&my_f,(void *) myparams, myparams->hessianh,m);

    //Invert the matrix
    gsl_permutation *p = gsl_permutation_alloc(myparams->sysdim);
    gsl_linalg_LU_decomp(m,p,&i);
    gsl_linalg_LU_invert(m,p,minv);
    gsl_permutation_free(p);

    /*Initialize the Vectors*/
    PROTECT(uconf = NEW_NUMERIC(myparams->sysdim));
    PROTECT(est = NEW_NUMERIC(myparams->sysdim));
    PROTECT(lconf = NEW_NUMERIC(myparams->sysdim));
    p_uconf=NUMERIC_POINTER(uconf);
    p_est=NUMERIC_POINTER(est);
    p_lconf=NUMERIC_POINTER(lconf);
    /*Calculate Confidence Intervals*/
    for (i = 0; i < myparams->sysdim; i++)
    {
       tmp=exp( gsl_vector_get(v,i)
		    + 1.959964 *sqrt(gsl_matrix_get(minv,i,i)));
       if(!myparams->method) tmp=tmp*
				  pow(myparams->nmols,myparams->order[i] -1);  
       p_uconf[i]=tmp;
    }

    if(myparams->method) 
      {
	for (i = 0; i < myparams->sysdim; i++){ 
	  p_est[i]= exp(gsl_vector_get(v,i));
	}
      }
    else
      {
	for (i = 0; i < myparams->sysdim; i++){ 
	  p_est[i]= exp(gsl_vector_get(v,i))*
	    pow(myparams->nmols,myparams->order[i]-1);
	}
      }

    for (i = 0; i < myparams->sysdim; i++)
    {
        tmp=exp( gsl_vector_get(v,i)
		    - 1.959964 *sqrt(gsl_matrix_get(minv,i,i)));
	if(!myparams->method) 
	  tmp=tmp*pow(myparams->nmols,myparams->order[i] -1);
	p_lconf[i]=tmp;

    }

    gsl_matrix_free(m);
    /*
    if(printh==1){
	printf("Covariance matrix:\n");
	for (i = 0; i < myparams->sysdim; i++)
	{
	    for (j = 0; j < myparams->sysdim; j++)
	    {
		printf("%f ",gsl_matrix_get(minv,i,j));
	    }
	    printf("\n");
	}
    }
    */
    gsl_matrix_free(minv);    
    /* Create List & atach the relevant vectors*/
    PROTECT(list_names = allocVector(STRSXP, 3));
    for(i = 0; i < 3; i++)
      SET_STRING_ELT(list_names, i,  mkChar(names[i]));
    PROTECT(list = allocVector(VECSXP, 3)); //Alocate list's elements
    SET_VECTOR_ELT(list, 0, uconf);
    SET_VECTOR_ELT(list, 1, est);
    SET_VECTOR_ELT(list, 2, lconf);
    setAttrib(list, R_NamesSymbol, list_names); //Set list's names
    UNPROTECT(5);
    return(list);
}


/**
 * Optimization procedure
 * 
 * Minimizes the log-likelihood function using the simplex method.
 * If the simplex has converged prints Confidence intervals of the mles.
 *
 * @param myparams A pointer to a spar_ data structure.
 */
int optim(spar **mypar)
{
    spar * myparams=*mypar;
    size_t iter = 0;
    int status,i;
    double size;

    const gsl_multimin_fminimizer_type *T = 
#ifdef SIMPLEX2
    gsl_multimin_fminimizer_nmsimplex2rand;
    //gsl_multimin_fminimizer_nmsimplex2;
#else
    gsl_multimin_fminimizer_nmsimplex;
#endif
    gsl_multimin_fminimizer *s = NULL;
    gsl_multimin_function minex_func;

            
    /* Initialize method and iterate */
    minex_func.n = myparams->sysdim;
    minex_func.f = &my_f;
    minex_func.params = (void *)myparams;


    s = gsl_multimin_fminimizer_alloc (T, myparams->sysdim);
    gsl_multimin_fminimizer_set (s, &minex_func, myparams->x, myparams->ss);

    do
    {
	iter++;
	status = gsl_multimin_fminimizer_iterate(s);

	if (status) 
	    break;

	size = gsl_multimin_fminimizer_size (s);
	status = gsl_multimin_test_size (size, myparams->simplexsize);

	if (status == GSL_SUCCESS)
	{
	    printf ("after %5d iterations converged to log-minimum at\n",iter);
	    for(i=0;i<myparams->sysdim;i++)
	    {
		printf("%6.6f ",exp(gsl_vector_get(s->x, i)));
	    }
	    printf("\n");
	    printf ("f() = %7.3f size = %1.3e\n", 
		    s->fval, size);
	}

	if(  iter % 100==0 )
	{
	    printf("%5d: ",iter);
	    for(i=0;i<myparams->sysdim;i++){
		printf("%3.3f ",exp(gsl_vector_get(s->x, i)));
	    }

	    printf ("f() = %7.5f log-size = %.3f\n", 
		    s->fval, size);
	}
    }
    while (status == GSL_CONTINUE && iter < myparams->maxiter);
    //turn parameters to exp scale
    for(i=0;i<myparams->sysdim;i++)
      {
	//gsl_vector_set(myparams->x,i, exp(gsl_vector_get(s->x, i)));
	gsl_vector_set(myparams->x,i, gsl_vector_get(s->x, i));
      }

    gsl_multimin_fminimizer_free(s);
    return 0;
}

double proflik(spar ** mypar)
{
    spar * myparams=*mypar;
    size_t iter = 0;
    int status,i,j;
    double size;
    const int sysdim=myparams->sysdim-1;

    gsl_vector *x = gsl_vector_alloc(sysdim);
    gsl_vector *ss = gsl_vector_alloc(sysdim);
    for(i=0;i<sysdim;i++)
    {
	
	if(myparams->pindex!=i)
	{ 
	    gsl_vector_set(x,i,gsl_vector_get(myparams->x,i));
	    gsl_vector_set(ss,i,gsl_vector_get(myparams->ss,i));
	} else
	{
	    for(j=i;j<sysdim;j++)
	    {
		gsl_vector_set(x,j,gsl_vector_get(myparams->x,j+1));
		gsl_vector_set(ss,j,gsl_vector_get(myparams->ss,j+1));
	    }
	    i=sysdim;
	}
    }
#ifdef SIMPLEX2
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
#else
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
#endif
    gsl_multimin_fminimizer *s = NULL;
    gsl_multimin_function minex_func;

            
    /* Initialize method and iterate */
    minex_func.n = sysdim;
    minex_func.f = &my_pl;
    minex_func.params = (void *)myparams;

    s = gsl_multimin_fminimizer_alloc (T, sysdim);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do
    {
	iter++;
	status = gsl_multimin_fminimizer_iterate(s);

	if (status) 
	    break;

	size = gsl_multimin_fminimizer_size (s);
	status = gsl_multimin_test_size (size, myparams->simplexsize);

	

	if(  iter % 50==0 && myparams->verbose==1)
	{
	    printf("%5d: ",iter);
	    for(i=0;i<x->size;i++){
		printf("%3.3f ",exp(gsl_vector_get(s->x, i)));
	    }

	    printf ("f() = %7.3f log-size = %.3f\n", 
		    s->fval, size);
	}
    }
    while (status == GSL_CONTINUE && iter < myparams->maxiter);
    if (status != GSL_SUCCESS)
    {
	fprintf(stderr,"for theta[%d] the proflik did not converge\n",
		myparams->pindex);
    }
    for(i=0;i<sysdim;i++) 
    {	
	if(myparams->pindex!=i)
	{ 
	    gsl_vector_set(myparams->x,i,gsl_vector_get(s->x,i));
	} else
	{
	    for(j=i;j<sysdim;j++)
	    {
		gsl_vector_set(myparams->x,j+1,gsl_vector_get(s->x,j));
	    }
	    i=sysdim;
	}
    }
    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(ss);
    return(s->fval);
}
double proflik2(spar ** mypar)
{
    spar * myparams=*mypar;
    size_t iter = 0;
    int status,i,j=0;
    double size;
    const int sysdim=myparams->sysdim-2;

    gsl_vector *x = gsl_vector_alloc(sysdim);
    gsl_vector *ss = gsl_vector_alloc(sysdim);
    for(i=0;i<sysdim;i++)
    {
	if(myparams->pindex==i || myparams->pindex2==i) j++;
	gsl_vector_set(x,i,gsl_vector_get(myparams->x,j));
	gsl_vector_set(ss,i,gsl_vector_get(myparams->ss,j));
	j++;
    }
#ifdef SIMPLEX2
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
#else
    const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex;
#endif
    gsl_multimin_fminimizer *s = NULL;
    gsl_multimin_function minex_func;

            
    /* Initialize method and iterate */
    minex_func.n = sysdim;
    minex_func.f = &my_pl2;
    minex_func.params = (void *)myparams;

    s = gsl_multimin_fminimizer_alloc (T, sysdim);
    gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

    do
    {
	iter++;
	status = gsl_multimin_fminimizer_iterate(s);

	if (status) 
	    break;

	size = gsl_multimin_fminimizer_size (s);
	status = gsl_multimin_test_size (size, myparams->simplexsize);

	

	if(  iter % 50==0 && myparams->verbose==1)
	{
	    printf("%5d: ",iter);
	    for(i=0;i<x->size;i++){
		printf("%3.3f ",exp(gsl_vector_get(s->x, i)));
	    }

	    printf ("f() = %7.3f log-size = %.3f\n", 
		    s->fval, size);
	}
    }
    while (status == GSL_CONTINUE && iter < myparams->maxiter);
    if (status != GSL_SUCCESS)
    {
	fprintf(stderr,"for theta[%d] the proflik did not converge\n",
		myparams->pindex);
    }
    /*
    for(i=0;i<sysdim;i++) 
    {	
	if(myparams->pindex!=i)
	{ 
	    gsl_vector_set(myparams->x,i,gsl_vector_get(s->x,i));
	} else
	{
	    for(j=i;j<sysdim;j++)
	    {
		gsl_vector_set(myparams->x,j+1,gsl_vector_get(s->x,j));
	    }
	    i=sysdim;
	}
    }
    */
    gsl_multimin_fminimizer_free(s);
    gsl_vector_free(x);
    gsl_vector_free(ss);
    return(s->fval);
}

/*
  Optimization based on BFGS
*/

int optim2(spar **mypar)
{
    spar * myparams=*mypar;
    size_t iter = 0;
    int status,i;
    
    const gsl_multimin_fdfminimizer_type *T;
    gsl_multimin_fdfminimizer *s;
    gsl_multimin_function_fdf minex_func;

            
    /* Initialize method and iterate */
    minex_func.n = myparams->sysdim;
    minex_func.f = &my_f;
    minex_func.df = &numgrad;
    minex_func.fdf= &my_fdf;
    minex_func.params = (void *)myparams;

    //T = gsl_multimin_fdfminimizer_conjugate_fr;
    //T = gsl_multimin_fdfminimizer_conjugate_pr;
    //T = gsl_multimin_fdfminimizer_vector_bfgs2;
    T = gsl_multimin_fdfminimizer_vector_bfgs;
    s = gsl_multimin_fdfminimizer_alloc (T, myparams->sysdim);
    gsl_multimin_fdfminimizer_set (s, &minex_func, myparams->x
				   ,0.1 /*Initial Step size*/ 
				   ,.01); /* Tolerance   */

    do
    {
	iter++;
	status = gsl_multimin_fdfminimizer_iterate(s);

	if (status)
	  {
	    printf("Exiting with status: %d\n",status);
	    printf ("f() = %7.3f\n", s->f);
	    break;
	  }
	status = gsl_multimin_test_gradient (s->gradient, myparams->simplexsize);

	if (status == GSL_SUCCESS)
	{
	    printf ("after %5d iterations converged to log-minimum at\n",iter);
	    for(i=0;i<myparams->sysdim;i++)
	    {
		printf("%6.6f ",exp(gsl_vector_get(s->x, i)));
	    }
	    printf("\n");
	    printf ("f() = %7.3f\n", 
		    s->f);
	    if(myparams->printci==1)
	    {
	      /*Print Covariance matrix*/
	      varcov(&myparams,s->x,0);
	    }
	    break;
	}

	if(  iter % 1==0 )
	{
	    printf("%5d: ",iter);
	    for(i=0;i<myparams->sysdim;i++){
		printf("%3.3f ",exp(gsl_vector_get(s->x, i)));
	    }

	    printf ("f() = %7.5f, ||df|| = %7.5f\n", s->f,
		    gsl_blas_dnrm2(s->gradient));
	}
    }
    while (status == GSL_CONTINUE && iter < myparams->maxiter);
    
    for(i=0;i<myparams->sysdim;i++)
      {
	gsl_vector_set(myparams->x,i, gsl_vector_get(s->x, i));
      }
        
    gsl_multimin_fdfminimizer_free(s);
    return 0;
}

#ifdef HAVEOOL
/**
 * Use the Open Optimization Library that provides:
 * SPG — This is the "Spectral Projected Gradient Method"
 * optimizes under constrains: exp(-5) < theta[i] < NMOLS
 * /
void iteration_echo( ool_conmin_minimizer *M )
{
   double f = M->f;
   size_t ii, nn;

   nn = 8;

   for( ii = 0; ii < nn; ii++ )
       printf( "%6.4f, ",exp( gsl_vector_get( M->x, ii )));
   printf( ": %3.3f , %2.3f\n", f,ool_conmin_minimizer_size( M ) );
}
int optim3(spar **mypar)
{
    spar * myparams=*mypar;
    size_t iter = 0;
    int status,i;
    const ool_conmin_minimizer_type *T = ool_conmin_minimizer_spg; 
    ool_conmin_spg_parameters P;
    ool_conmin_function   F;
    ool_conmin_constraint C;
    ool_conmin_minimizer *M;
            
    /* Initialize method and iterate */
    F.n   = myparams->sysdim;
    F.f   = &my_f;
    F.df  = &numgrad;
    F.fdf = &my_fdf;
    F.params = (void *)myparams;

    C.n = myparams->sysdim;
    C.L = gsl_vector_alloc( C.n );
    C.U = gsl_vector_alloc( C.n );
    gsl_vector_set_all( C.L, -5.0 );
    gsl_vector_set_all( C.U,  log(myparams->nmols) );

    M = ool_conmin_minimizer_alloc( T, myparams->sysdim );
    ool_conmin_parameters_default( T, (void*)(&P) );
    ool_conmin_minimizer_set( M, &F, &C, myparams->x, (void*)(&P) );
    iteration_echo ( M );
    status=OOL_CONTINUE;
    do
    {
	iter++;
	ool_conmin_minimizer_iterate( M );
	status = ool_conmin_is_optimal( M );

	if (status != OOL_CONTINUE )
	{
	    printf ("after %5d iterations converged to log-minimum at\n",iter);
	    for(i=0;i<myparams->sysdim;i++)
	    {
		printf( "%6.9f, ",exp( gsl_vector_get( M->x, i )));

	    }
	    break;
	}

	if(  iter % 1==0 )
	{
		iteration_echo( M );
	}
    }
    while (status == OOL_CONTINUE && iter < myparams->maxiter);
        
    gsl_vector_free( C.L );
    gsl_vector_free( C.U );
    for(i=0;i<myparams->sysdim;i++) 
	    gsl_vector_set(myparams->x,i,exp(gsl_vector_get(M->x,i)));
    ool_conmin_minimizer_free( M );
    return 0;
}
#endif
