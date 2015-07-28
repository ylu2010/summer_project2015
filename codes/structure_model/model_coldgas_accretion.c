/*
 * model_coldgas_accretion.c
 *
 *  Created on: Jun 14, 2012
 *      Author: luyu
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_rng.h>


#include "variables.h"
#include "proto.h"
#include "cosmo.h"


void cold_gas_accretion(struct galaxy *gal, double thubble, double dt)
{
	int i, nbin;
	double cr, bar, z, r, sig0;

	z = Redshift;

	//bar = gal->RateHaloAccretion * dt; // no real cooling calculation; the cooling rate just follows the halo accretion rate
	cr = cooling_rate(gal, thubble, dt, 1); // based on cooling calculation
	bar = cr * dt;

	gal->RateCooling = cr;
	// masses
	gal->MassCold += bar;
	gal->MassHot -= bar;

	// distribution
	sig0 = gal->MassCold/(2.*M_PI * gal->RadiusDisc * gal->RadiusDisc);
	nbin = gal->nbin;
	for(i=0; i<nbin; i++)
	{
		r = 0.5*(gal->RadiusInner[i] + gal->RadiusOuter[i]);
		gal->SDensityCold[i] = sig0 * exp(-r/gal->RadiusDisc);
		//printf("%d %g %g %g %g %g\n", i, gal->RadiusInner[i], gal->RadiusOuter[i], gal->MassProfHalo[i], gal->SDensityCold[i], gal->SDensityStar[i] );
	}

}

void cold_gas_accretion_surface(struct galaxy *gal, double thubble, double dt)
{
	int i, j, k, nbin;
	double cr, bar, z, r, sig0, tdyn, t0, hr=0.0;
	double zhot;

	if (gal->MassHot > 0.0) zhot = gal->MetalHot/gal->MassHot;
	else zhot = 0.0;

	z = Redshift;
	tdyn = gal->RadiusHalo*UnitLength_in_cm*1e-5/gal->VelocityVirial/UnitTime_in_Second;

	//cr = 0.17*gal->RateHaloAccretion;  // no real cooling calculation; the cooling rate just follows the halo accretion rate
	//if(z<2.01) cr = 0.0 * cr;//0.17 * Par.BaryonAccretionFraction * gal->RateHaloAccretion * do_preheating(gal, z);
	cr = cooling_rate_shell(gal, thubble, dt, 1);  // based on cooling calculation

	//hr = 5.0005*gal->RateStarFormation*(200./gal->VelocityVirial)*(200./gal->VelocityVirial);
	//hr = 7e-5*gal->MassStar * 130*130;// * gal->VelocityVirial * gal->VelocityVirial;
	/*
	 * AGB heating Conroy et al. 2014 It does not produce a decreasing SFR as a function of time.
	 * The reason is it has a time dependence 1/t^1.25. Even with a weaker dependence 1/t, the predicted SFR is still too flat.
	 *
	j = floor(thubble/Bin_size_time);
	if (j > N_TIME_BIN) j=N_TIME_BIN;
	for(k=0, t0=0, hr=0; k<j; k++)
	{
		t0 = (thubble - (0.5+k) * Bin_size_time) * UnitTime_in_Megayears * 1e6;
		if (t0 > 1e7)
		{
			hr += gal->StarFormationHistory[k] / pow(t0, 1.25);
		}
		//printf("k=%d th=%g t0=%g sfr=%g\n", k, thubble, t0, gal->StarFormationHistory[k]);
	}
	hr =  0.5 * 8. * hr * UnitTime_in_Megayears * 1e6 * pow(gal->VelocityVirial/30.,2);
	printf("t=%g vel:%g j=%d hr=%g mstar=%g\n", thubble, gal->VelocityVirial, j, hr, gal->MassStar);
	*/
	hr = Par.GalaxyHeatingEfficiency * gal->MassStar * pow(gal->MassStar/2e10, 1.);// * pow(gal->VelocityVirial/120.,2);// * gal->MassStar/2e10;

	cr = dmax(cr - hr, 0.0);

	//if(z<2.01) cr = 0.0;
	//cr = 0.0;
	//cr += 0.5*gal->MassCloud / tdyn;
	gal->MassCloud -= gal->MassCloud / tdyn * dt;
	bar = cr * dt;

	gal->RateCooling = cr;
	// masses
	gal->MassCold += bar;
	gal->MassHot -= bar;

	//distribution
	sig0 = bar/(2.*M_PI * gal->RadiusDisc * gal->RadiusDisc);
	nbin = gal->nbin;
	for(i=0; i<nbin; i++)
	{
		r = 0.5*(gal->RadiusInner[i] + gal->RadiusOuter[i]);
		gal->SDensityCold[i] += sig0 * exp(-r/gal->RadiusDisc);
		//printf("%d %g %g %g %g %g %g\n", i, bar, gal->RadiusInner[i], gal->RadiusOuter[i], gal->MassProfHalo[i], gal->SDensityCold[i], gal->SDensityStar[i] );
        gal->SDensityCAR[i] = sig0 * exp(-r/gal->RadiusDisc);
		gal->MassMetalCold[i] += sig0 * exp(-r/gal->RadiusDisc) * zhot * M_PI * (gal->RadiusOuter[i]*gal->RadiusOuter[i] - gal->RadiusInner[i] * gal->RadiusInner[i]);
	}

	// metallicity

	gal->MetalCold += cr * dt * zhot;
	gal->MetalHot -= cr * dt * zhot;
}
