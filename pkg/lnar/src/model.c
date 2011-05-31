#include <math.h>
#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <string.h>
#include <gsl/gsl_blas.h>
#include <getopt.h>
#include "model.h"
#include "optim.h"
#include "integrate.h"

/*void derivs2(double t, const double *y, double *f, void *params)
{
    derivs(t,y,f,params);
}
*/
SEXP runmodel(SEXP ptrf, SEXP NPARAMS, SEXP NODES, SEXP maxiter,
	SEXP nmols, SEXP tcrit, SEXP hessianh, SEXP relerr,
	SEXP abserr, SEXP dataset, SEXP nr,SEXP nc, SEXP ord,
        SEXP thetas, SEXP usebfgs, SEXP method){

    SEXP list;
    int ny=LENGTH(dataset);
    int nth=LENGTH(thetas);
    int nor=LENGTH(ord);
    /*
    //Debug mesages
    Rprintf("size=%d\n",ny);
    Rprintf("nr=%d\n",INTEGER(nr)[0]);
    Rprintf("nc=%d\n",INTEGER(nc)[0]);
    */
    double *ytmp = (double *) R_alloc( ny, sizeof(double));
    double *thtmp = (double *) R_alloc( nth, sizeof(double));
    int *ortmp = (int *) R_alloc( nor, sizeof(int));
    double * p_thetas = NUMERIC_POINTER(thetas);
    double logbarp=0,vartmp;
    char *parfile=NULL,*outfile=NULL,*partruef=NULL;
    //process options
    int i,j;
    int printci=1; /* Defaults to not printing CI */
    int bfgs;

    
    for (j = 0; j < ny; j++) ytmp[j] = REAL(dataset)[j];
    for (j = 0; j < nth; j++) thtmp[j] = REAL(thetas)[j];
    for (j = 0; j < nth; j++) ortmp[j] = INTEGER(ord)[j];


    

    spar * myparams = (spar *) R_alloc(1,sizeof(spar));
    /* define constants  */
    myparams->dataset=ytmp;
    myparams->thetas=p_thetas;
    myparams->nrows=INTEGER(nr)[0];
    myparams->ncols=INTEGER(nc)[0];
    myparams->nmols=INTEGER(nmols)[0];
    myparams->sysdim=INTEGER(NPARAMS)[0];
    myparams->maxiter=INTEGER(maxiter)[0];
    myparams->simplexsize=REAL(tcrit)[0];
    myparams->hessianh=REAL(hessianh)[0];
    myparams->relerr=REAL(relerr)[0];
    myparams->abserr=REAL(abserr)[0];
    myparams->odesize=INTEGER(NODES)[0];
    myparams->thetatrue=NULL;
    myparams->f= (funderivs ) R_ExternalPtrAddr(ptrf);

    myparams->jac=NULL;
    myparams->printci=printci;
    myparams->order = ortmp;
    myparams->method=INTEGER(method)[0];
    bfgs=INTEGER(usebfgs)[0];


    switch(myparams->method) {
	case 0:
	    Rprintf("Restarting Method - Molecules\n");
	    break;
        case 1:
	    Rprintf("Restarting Method - Concentrations\n");
	    break;
	case 4: 
	    Rprintf("Non-Restarting Method\n");
	    break;
	default:
	    Rprintf("Please specify a valid method\n");
	    exit;
    }

    /* Parameters Starting point */
    myparams->x = gsl_vector_calloc (myparams->sysdim);
    for(i=0;i<myparams->sysdim;i++)
    {
	gsl_vector_set(myparams->x,i,log(thtmp[i]));
    }

    /* Set initial step sizes to 1: correspond to exp(1) */
    myparams->ss = gsl_vector_alloc (myparams->sysdim);
    gsl_vector_set_all (myparams->ss, 1);
    /* Define initial simplex size or convergence criterion */
    if(bfgs==1 ) {
	myparams->simplexsize= 
	  ( REAL(tcrit)[0]>0 ?  REAL(tcrit)[0]: 1) ;
    }
    else {
    if( tcrit>0) gsl_vector_set_all(myparams->ss, REAL(tcrit)[0]);
    }    


    if(myparams->method==0)
      {
	for (i=0; i<myparams->sysdim; i++)
	  {
	    gsl_vector_set(myparams->x,i,
			   gsl_vector_get(myparams->x,i)+ 
			   (1 - myparams->order[i])*log(myparams->nmols)
			   );
	  }
      }
    

    myparams->logc=0;   
    if(bfgs) optim2(&myparams); else  optim(&myparams);
    list=varcov(&myparams,myparams->x,0);
    gsl_vector_free(myparams->x);
    gsl_vector_free(myparams->ss);
    //return R_NilValue;
    return(list);
}

SEXP calclik(SEXP ptrf, SEXP NPARAMS, SEXP NODES,
	SEXP nmols, SEXP relerr,SEXP abserr, 
	SEXP dataset, SEXP nr,SEXP nc, SEXP ord,
        SEXP thetas, SEXP method){

    SEXP ans;
    int ny=INTEGER(NODES)[0];
    int nth=LENGTH(thetas);
    int nor=LENGTH(ord);
    double * p_thetas = NUMERIC_POINTER(thetas);
    double * p_dataset= NUMERIC_POINTER(dataset);
    double logbarp=0,vartmp;
    int * p_order = INTEGER_POINTER(ord);
    char *parfile=NULL,*outfile=NULL,*partruef=NULL;
    //process options
    int i,j;
    

    spar * myparams = (spar *) R_alloc(1,sizeof(spar));
    /* define constants  */
    myparams->dataset=p_dataset;
    myparams->thetas=p_thetas;
    myparams->nrows=INTEGER(nr)[0];
    myparams->ncols=INTEGER(nc)[0];
    myparams->nmols=INTEGER(nmols)[0];
    myparams->sysdim=INTEGER(NPARAMS)[0];
    myparams->relerr=REAL(relerr)[0];
    myparams->abserr=REAL(abserr)[0];
    myparams->odesize=INTEGER(NODES)[0];
    myparams->thetatrue=NULL;
    myparams->f= (funderivs ) R_ExternalPtrAddr(ptrf);

    myparams->jac=NULL;
    myparams->order = p_order;
    myparams->method=INTEGER(method)[0];
    
    switch(myparams->method) {
	case 0:
	    Rprintf("Restarting Method - Molecules\n");
	    break;
        case 1:
	    Rprintf("Restarting Method - Concentrations\n");
	    break;
	case 4: 
	    Rprintf("Non-Restarting Method\n");
	    break;
	default:
	    Rprintf("Please specify a valid method\n");
	    exit;
    }

    myparams->logc=0;
    PROTECT(ans = NEW_NUMERIC(1));
    REAL(ans)[0]=llik(&myparams, p_thetas);    
    UNPROTECT(1);
    //return R_NilValue;
    return(ans);
}
