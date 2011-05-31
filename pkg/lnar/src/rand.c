/* $Id: rand.c,v 1.1 2009/02/01 20:05:32 linus Exp linus $
 * Generates a random seed for gsl library
 * */

#include "rand.h"

unsigned int arandom_seed()
{
    unsigned int seed;
    struct timeval tv;
    FILE *devrandom;

    if ((devrandom = fopen("/dev/random","r")) == NULL) {
        gettimeofday(&tv,0);
        seed = tv.tv_sec + tv.tv_usec;
        #ifdef DEBUG
            printf("Got seed %u from gettimeofday()\n",seed);
        #endif
    } else {
        fread(&seed,sizeof(seed),1,devrandom);
        #ifdef DEBUG
            printf("Got seed %u from /dev/random\n",seed);
        #endif
        fclose(devrandom);
    }    
    return(seed);
}

void rand_init(unsigned int * seed){
	/* initialize the rng from enviromentals */
	gsl_rng_env_setup(); 
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	if(seed==NULL) {
		/*Peeks a seed pseudo-randomly */
		gsl_rng_set (r, arandom_seed());
	}		/*End of the initialization */
	else{
		gsl_rng_set(r, *seed);
	}
}
void rand_end(){
	gsl_rng_free(r);
}
