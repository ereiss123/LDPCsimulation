/* rand_gsl.h - Random number generators. */
/* Formatted to be compatible with Neal's 
   rand.h but using GSL random generators. 

   This version:
   By Chris Winstead, Utah State University.

   Original rand.h copyright message:
   Copyright (c) 1992 by Radford M. Neal 
*/

// GSL includes:                        
#include <gsl/gsl_rng.h>                          
#include <gsl/gsl_randist.h>                      


/* SET RANDOM NUMBER SEED. */

gsl_rng * rng_engine;

#define ran_seed(s) rng_engine = gsl_rng_alloc (gsl_rng_taus);	\
  gsl_rng_set (rng_engine, s)


/* GENERATE RANDOM NUMBERS. */

#define ranf() \
  gsl_rng_uniform(rng_engine)

#define ranu() \
  gsl_rng_uniform_pos(rng_engine);

#define rani(n) \
  gsl_rng_uniform_int(rng_engine,n)

#define rann() \
  gsl_ran_ugaussian(rng_engine)

#define rane() \
  (-log(ranu()))		                  /* From exponential */

#define ranc() \
  (tan(3.141592654*(ranu()-0.5)))		                      /* From Cauchy */
