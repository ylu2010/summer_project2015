/*
 * main_run_pre_sample.c
 *
 *  Created on: May 29, 2014
 *      Author: luyu
 *     main function run a parameter space sampler to produce parameter samples and run the model using the parameter samples.
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "variables.h"
#include "proto.h"
#include "ihs.hpp"

int set_varying_parameters(double *params, int nparams)
{
	if(Do_preheating)
	{
		Par.PreheatEntropy = params[0];
		Par.PreheatEntropySlope = params[1];
	}
	else
	{
		Par.SNLoadingFactor = params[0];
		Par.SNLoadingFactorIndex = params[1];
	}
    //Par.GalaxyHeatingEfficiency = params[2];
    //Par.Yield = params[3];
    //Par.ZFractionYieldToHot = params[4];
	Par.DiskRadiusFactor = params[2];

	return 0;
}

int set_predictions(double *p, int np, struct galaxy *gal, int ihalo)
{
	/*
    double mdust_pee, mzcold;
    p[0+ihalo*5] = log10(gal->MassHalo);
    p[1+ihalo*5] = log10(gal->MassStar);
    p[2+ihalo*5] = log10(gal->MassColdAtomic/gal->MassStar);
    mdust_pee = pow(10.,(0.86*log10(gal->MassStar)-1.31));
    mzcold = dmax(0.0, gal->MetalCold - mdust_pee);
    p[3+ihalo*5] = log10(gal->MetalCold/gal->MassCold);
    p[4+ihalo*5] = log10(gal->RadiusHalfStar)+3.; //change unit to kpc from mpc
    */
    return 0;
}

int general_parameter_sample(int nparams, int point_num, struct interval *param_interval, int seed, double *p)
{
	int ipoint, iparam;
	int *x;

	x = (int *) malloc(nparams * point_num * sizeof(int));

	sample_parameters(x, nparams, point_num, seed);

	for(ipoint=0; ipoint<point_num; ipoint++)
	{
		for(iparam=0; iparam<nparams; iparam++)
	    {
	    	p[ipoint*nparams+iparam] = (((double) x[ipoint*nparams+iparam] - 0.5) / point_num ) *
	    			(param_interval[iparam].max-param_interval[iparam].min) + param_interval[iparam].min;
	    }
	}

	free(x);
	return 0;
}


int main( int argc, const char* argv[] )
{
	int nparams = 3;
	int npreds = 5*2;
	int point_num = 300;
	int seed = 17;
	int ihalo = 0;
	int mode, write_pred_saparately;
	double dbuf;
	double *params, *preds;
	double *p;
	int irun, iparam;
	struct interval *param_interval;
	char input_param_list[200];
	FILE *fp;

    if(argc < 3)
    {
        printf("\n  usage: sample_gstructure 0 <random_seed (cannot be 0!)>\n\n");
        printf("\n  usage: sample_gstructure 1 <ParameterList>\n\n");
        exit(1);
    }

    mode = atoi(argv[1]);

    param_interval = (struct interval *) malloc(sizeof(struct interval)*nparams);
    if (mode == 0)
    {
    	seed = atoi(argv[2]);
    	printf("In the sampling mode... The code is generating parameter samples with a random seed=%d...\n\n", seed);

    	param_interval[0].min = 0.0;
    	param_interval[0].max = 10.0;
    	param_interval[1].min = 0.0;
    	param_interval[1].max = 4.0;
    	param_interval[2].min = 0.0;
    	param_interval[2].max = 1.0;
    	if(!(  fp=fopen("param_range.txt","w")))
    	{
    		printf("I can not open the file param_range.txt.");
    		exit(1);
    	}
    	fprintf(fp, "2 %d\n", nparams);
    	for(iparam=0; iparam<nparams; iparam++)
    	{
    		fprintf(fp, "%g %g\n", param_interval[iparam].min, param_interval[iparam].max);
    	}
    	fclose(fp);

    	p = (double *) malloc(nparams * point_num * sizeof(double));
    	general_parameter_sample(nparams, point_num, param_interval, seed, p);
    }
    else
    {
    	if(argc < 3)
    	{
    		printf("In a mode where a pre-generated parameter list is expected.\n");
    		printf("\n  usage: sample_gstructure 1 <ParameterList>\n\n");
    	    exit(1);
    	}
    	strcpy(input_param_list, argv[2]);
    	if(!(  fp=fopen(input_param_list,"r")))
    	{
    	    printf("I cannot open the file %s.", input_param_list);
    	    exit(1);
    	}
    	printf("Reading file %s.", input_param_list);
    	fscanf(fp, "%d %d\n", &point_num, &nparams);
    	printf("%d %d\n", point_num, nparams);

    	p = (double *) malloc(nparams * point_num * sizeof(double));
    	for(irun=0; irun<point_num; irun++)
    	{
    		for(iparam=0; iparam<nparams; iparam++)
    		{
    			fscanf(fp, "%lf ", &dbuf);
    			p[irun*nparams+iparam] = dbuf;
    			//printf("%g %g ", dbuf, p[irun*nparams+iparam]);
    		}
    		//printf("\n");
    	}
    	fclose(fp);

    }
	params = (double *) malloc(sizeof(double) * nparams);
	preds = (double *) malloc(sizeof(double) * npreds);

    setup_run();

	fp_list=fopen("list.dat","w");

    for(irun=0; irun<point_num; irun++)
    {
		fprintf(fp_list, "%d ", irun);
    	for(iparam=0; iparam<nparams; iparam++)
    	{
    		params[iparam] = p[irun*nparams+iparam];
			fprintf(fp_list, "%f ", params[iparam]);
    	}
		fprintf(fp_list, "\n");

    	run_galaxy(params, nparams, preds, npreds, mode=0, irun);
    }
    finalize_run();

    free(preds);
    free(params);
    free(p);
    free(param_interval);
	fclose(fp_list);

    return 0;
}


