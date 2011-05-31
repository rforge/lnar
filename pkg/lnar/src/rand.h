/* $Id: rand.h,v 1.2 2009/02/02 02:29:24 linus Exp linus $
 * Generates a random seed for gsl library
 * */

#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdlib.h>

gsl_rng * r; 
unsigned int arandom_seed();
const gsl_rng_type * T;
void rand_init(unsigned int * seed);
void rand_end();
