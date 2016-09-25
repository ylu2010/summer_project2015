/*
 * misc.c
 *
 *  Created on: Mar 19, 2013
 *      Author: luyu
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>


#include "variables.h"
#include "proto.h"


#define MAXBLOCKS 256

double dmax(double x, double y)
{
    if(x > y)
        return x;
    else
        return y;
}


double dmin(double x, double y)
{
    if(x < y)
        return x;
    else
        return y;
}


/* returns the number of cpu-ticks in seconds that
 * have elapsed. (or the wall-clock time)
 */
double second(void)
{
  return ((double)((unsigned int)clock()))/CLOCKS_PER_SEC;

  /* note: on AIX and presumably many other 32bit systems,
   * clock() has only a resolution of 10ms=0.01sec
   */
}
/* returns the time difference between two measurements
 * obtained with second(). The routine takes care of the
 * possible overflow of the tick counter on 32bit systems.
 */
double timediff(double t0,double t1)
{
  double dt;

  dt=t1-t0;

  if(dt<0)  /* overflow has occured */
    {
      dt=t1 + pow(2,32)/CLOCKS_PER_SEC - t0;
    }

  return dt;
}


const gsl_rng_type *rgtype, *rutype;
gsl_rng *rg, *ru;

void init_random_number_gstructure(void)
{
    gsl_rng_env_setup();
    gsl_rng_default_seed = 3;

    rgtype = gsl_rng_default;
    rg = gsl_rng_alloc (rgtype);

    rutype = gsl_rng_default;
    ru = gsl_rng_alloc (rutype);
}

void free_random_number_gstructure(void)
{
    gsl_rng_free(rg);
    gsl_rng_free(ru);
}

double random_number_gstructure(int flag)
{
    double rn;

    if(flag == 0) rn = gsl_rng_uniform(ru);

    if(flag == 1) rn = gsl_ran_ugaussian(rg);
//printf("flag=%d rn=%g\n", flag, rn);
    return rn;
}

double interpolate_bipoint(double *x, double *y, int n, double x0, int flag)
{
    int i;
    double x1, x2, y0;


    if (flag >=0 )
    {
    	x1 = x[0]; x2 = x[1];
    	i=0;
    	do
    	{
    		i++;
    		//printf("interpolate: i=%d %g %g\n", i, x[i], x0);
    		if(i>n-1)
    		{
    			fprintf(stderr,"the interpolation point is out of range! i=%d x1=%g x2=%g x?=%g\n", i, x1, x2, x0);
    			exit(0);
    		}
    	}while(x[i] < x0);
    }
    else
    {
    	x1 = x[n-1]; x2 = x[n-2];
    	i=n-1;
    	do
    	    {
    	    	i--;
    	    	//printf("interpolate: i=%d %g %g\n", i, x[i], x0);
    	    	if(i<0)
    	    	{
    	    		fprintf(stderr,"the interpolation point is out of range! i=%d x1=%g x2=%g x?=%g\n", i, x1, x2, x0);
    	    		exit(0);
    	    	}
    	    }while(x[i] < x0);
    }
    if(x[i+1] != x[i])
    	y0 = (y[i+1] - y[i])/(x[i+1] - x[i]) * (x0 - x[i]) + y[i];
    else
    	y0 = 0.5* (y[i]+ y[i+1]);
    //printf("int: n=%d i=%d %g %g\n", n, i, x[i], x[i+1]);
    //fprintf(stderr,"the interpolation point is out of range! %d %d x1=%g x2=%g x?=%g %g %g y?=%g\n", flag, i , x[i], x[i+1], x0, y[i], y[i+1], y0);
    return y0;
}

