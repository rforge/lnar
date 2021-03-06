#include "integrate.h"
#include <stdio.h>
#include <math.h>
//#include <malloc.h>
#include <stdlib.h>
#include <string.h>

/* Handle errors within R */
#include <R.h>

extern void integrate(double *y,double *t,double * tout,
    double * theta,int  * odesize,
    void (*devis)(), double * reltol, double * abstol)
{
  int lwr= (*odesize) * (*odesize) +22 + (*odesize) * 15;
  int liw =(*odesize)+20;
  double          rwork[lwr];
  double  atol,  rtol;
  /*double          rtol[odesizep],y_[odesizep]; */
  int             iwork[liw];
  int itol = 1; //Scalar
  int  itask, istate, iopt, jt, i;

  iwork[5]=8000;


  rtol= *reltol;
  atol= *abstol;

  itask = 1;
  istate = 1;
  iopt = 0;
  jt = 2;
  lsoda_(devis, odesize, y, t, tout, &itol, &rtol, &atol, 
      &itask, &istate, &iopt, rwork, &lwr, iwork, &liw, NULL, 
      &jt,  theta);


  if (istate <= 0) {
    error("error istate = %d\n", istate);
  }
}
