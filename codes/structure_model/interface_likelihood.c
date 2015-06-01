/*
 * interface_likelihood.c
 *
 *  Created on: May 29, 2014
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
    Par.GalaxyHeatingEfficiency = params[2];
    Par.Yield = params[3];
    Par.ZFractionYieldToHot = params[4];

	return 0;
}

int set_predictions(double *p, int np, struct galaxy *gal, int ihalo)
{
    double mdust_pee, mzcold;
    p[0+ihalo*5] = log10(gal->MassHalo);
    p[1+ihalo*5] = log10(gal->MassStar);
    //p[2+ihalo*5] = log10(gal->MassColdAtomic/gal->MassStar);
    p[2+ihalo*5] = log10(gal->MassCold/gal->MassStar);
    mdust_pee = pow(10.,(0.86*log10(gal->MassStar)-1.31));
    mzcold = dmax(0.0, gal->MetalCold - mdust_pee);
    p[3+ihalo*5] = log10(gal->MetalCold/gal->MassCold);
    p[4+ihalo*5] = log10(gal->RadiusHalfStar)+3.; //change unit to kpc from mpc
    return 0;
}


void likelihood_init(void)
{
	setup_run();
}

void likelihood_finalize(void)
{
	finalize_run();
}

double chi2_mstar_mhalo(double lgms, double lgmh)
{
	double lgfs;
	double val;
	lgfs = lgms - lgmh;
	if(lgmh > 11.5) val=(lgfs+1.58184)/(0.269668*0.5);
	else val=(lgfs+2.34598)/(0.289756*0.5);
	return val*val;
}

double chi2_mstar_mhalo_z2(double lgms, double lgmh)
{
        double lgfs;
        double val;
        lgfs = lgms - lgmh;
        if(lgmh > 11.5) val=(lgfs+1.68323)/(0.269668*0.5);
        else val=(lgfs+2.70315)/(0.289756*0.5);
        return val*val;
}


double chi2_fcold_mstar(double lgfc, double lgms)
{
	double val;
	double lgfgas, sig_fgas;

	lgfgas=-0.48*lgms+4.39; //Peeples et al. 2014
	//lgfgas = -0.43*lgms+3.89; // Papastergis et al. 2012 (from Peeples et al. 2014)
	sig_fgas=0.5;
	val=(lgfc-lgfgas)/sig_fgas;

	return val*val;
}

double chi2_zcold_mstar(double lgzc, double lgms)
{
	double val;
	double lgnoxy, lgzcold_maiolino;
	double lgzcold_obs, sig_zcold;

	double lgm0=11.18;
	double k0=9.04;
	double a0=-0.0864;
	lgnoxy=a0*pow(lgms-lgm0,2)+k0-12.;
	lgzcold_maiolino=lgnoxy + log10(15.999/1.366/0.44);

	lgzcold_obs=lgzcold_maiolino;
	sig_zcold=0.5;
	val = (lgzc-lgzcold_obs)/sig_zcold;

	return val*val;
}

double chi2_rd_mstar(double lgrd, double lgms)
{
	double a, b, c;
	double lgrd_z0_dutton11, sig;
	double val;

	a = 0.18; b = 0.52; c=1.8;
	lgrd_z0_dutton11 = 0.72 + a*(lgms-10.44) + ((b-a)/c) * log10(0.5+0.5*pow(10,(lgms-10.44)*c));
	sig = 0.25;
	val = (lgrd - lgrd_z0_dutton11)/sig;

	return val*val;
}

double likelihood( int ndim, double *pos, int flag)
{
	int mode, itag, ihalo, nhalo=2;
	int npreds = 5*2;
	double preds[5*2];
	double lgmhalo, lgmstar, lgfcold, lgzcold, lgrd;

	double val=0.0;
	run_galaxy(pos, ndim, preds, npreds, mode=0, itag=0);
	//printf("lgmhalo=%g lgmstar=%g lgfcold=%g lgzcold=%g lgrd=%g \n", preds[0], preds[1], preds[2], preds[3], preds[4]);
	for(ihalo=0; ihalo<nhalo; ihalo++)
	{
		lgmhalo = preds[ihalo*5+0];
		lgmstar = preds[ihalo*5+1];
		lgfcold = preds[ihalo*5+2];
		lgzcold = preds[ihalo*5+3];
		lgrd = preds[ihalo*5+4];

		val += chi2_mstar_mhalo(lgmstar, lgmhalo);
		//val += chi2_mstar_mhalo_z2(lgmstar, lgmhalo);
		//printf("val1=%g\t", -0.5*val);
		val += chi2_fcold_mstar(lgfcold, lgmstar);
		//printf("val2=%g\t", -0.5*val);
		val += chi2_zcold_mstar(lgzcold, lgmstar);
		//printf("val3=%g\t", -0.5*val);
		//val += chi2_rd_mstar(lgrd, lgmstar);
		//printf("val4=%g\n", -0.5*val);
	}
	return -0.5*val;
}

int main( int argc, const char* argv[] )
{
	int ndim = 5;
	//double pos[4] = {14.411221 ,       0.1626767,      0.017949153 ,    0.2034860014};
	//double pos[4] = {29.388885   ,    0.10302432  ,    0.021399688   ,    0.65421559};
	//double pos[4] = {29.737126   ,    0.75297192  ,    0.034561284   ,    0.93016779};
	//double pos[4] = {14.259926   ,    0.41062301   ,   0.022712228    ,  0.081299162};
	//double pos[4] = {14.163654   ,    0.67518377    ,   0.02629324   ,    0.41984442};

	//double pos[4] = {0.97832452  ,     0.97671701   ,   0.038999137   ,    0.38295562};
	double pos[5] ={11.597175,       0.17512586,        1.0331305,      0.015474346,       0.72381624};

	double val;
	likelihood_init();
	val = likelihood( ndim, pos, 0);
	likelihood_finalize();
	printf("main: val=%g\n", val);
	return 0;
}

