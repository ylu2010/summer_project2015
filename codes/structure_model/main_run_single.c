/*
 * main_run_single.c
 *
 *  Created on: May 21, 2014
 *      Author: luyu
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables.h"
#include "proto.h"
//#include "ihs.hpp"

int set_varying_parameters(double *params, int nparams)
{
	return 0;
}

int set_predictions(double *p, int np, struct galaxy *gal, int ihalo)
{

	return 0;
}


int main( int argc, const char* argv[] )
{
	int nparams = 1;
	int npreds = 1;

	int ihalo;
	int mode;
	double *params, *preds;

    setup_run();

    run_galaxy(params, nparams, preds, npreds, 1, 0);

    finalize_run();

    return 0;
}
