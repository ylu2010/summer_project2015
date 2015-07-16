/*
 * run_galaxy.c
 *
 *  Created on: May 21, 2014
 *      Author: luyu
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "variables.h"
#include "proto.h"
#include "cosmo.h"

double UnitMass_in_g;
double UnitLength_in_cm;
double UnitTime_in_Megayears;
double UnitTime_in_Second;
double UnitVelocity_in_cm_per_s;
double BaryonFrac = 0.17;
//double SolarMetallicity = 0.0142;
double SolarMetallicity = 0.0134;
double MinimumMetallicityRelativeToSolar = 0.001;

int Metal_gas_evolu=1;
int Mah_simu = 0;
int Do_preheating = 0;
int Do_reinfall = 0;
int N_halo = 11;
float Mass_bin = 1.0;
float LogHaloMassArray[11]={10.0,10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0,12.25,12.5};
int Write_pred_file=1;
int Write_pred_saparately=0;
int Write_hist_file=1;
int Write_prof_file=1;
double Redshift;
double Redshift_end=0.0;

struct parameter Par;

FILE *fp_hist;
FILE *fp_pred;
FILE *fp_disc;
FILE *fp_list;

int setup_run(void)
{

	UnitMass_in_g = MSUN_IN_G;
	UnitLength_in_cm = PC_IN_CM*1e6;
	UnitTime_in_Megayears = 1e3 ;
	UnitTime_in_Second = 1e9*365*24*3600.;
	UnitVelocity_in_cm_per_s = 1e5;

	Par.Reionization_z0 = 12;
	Par.Reionization_zr = 11;

	Par.BaryonAccretionFraction = 1.0;
	Par.StellarMassLossFraction = 0.43;//0.35/1.35;
	Par.Yield = 0.03;

	Par.MassFractionEjectToHot = 0.0;

	if ( Do_preheating)
	{
		Par.PreheatEntropySlope = 0.2;
		Par.PreheatEntropy = 10;//17;//12;//* pow(pow(10.,Mass_bin)/1e12, Par.PreheatEntropySlope);
		Par.EntropyProfileIndex = 0.0;
		Par.DiskRadiusFactor = 0.7;
		Par.ZFractionYieldToEject = 0.0;// for PR model
		Par.ZFractionYieldToHot = 0.1;  // for PR model
		Par.GalaxyHeatingEfficiency = 1.;//1.1;
		Par.SNLoadingFactor= 1;
		Par.SNLoadingFactorIndex = 0.;
	}
	else
	{
		Par.PreheatEntropySlope = 0.;
		Par.PreheatEntropy = 0.0;
		Par.EntropyProfileIndex = 1.1;//4./3;
		Par.DiskRadiusFactor = 0.8;    // for EJ model
		Par.ZFractionYieldToEject = 0.;// for EJ model
		Par.ZFractionYieldToHot = 0.1;  // for EJ model
		Par.GalaxyHeatingEfficiency = 0.;
		Par.SNLoadingFactor= 1;
		Par.SNLoadingFactorIndex = 2;
		if (Do_reinfall)
		{
			Par.DiskRadiusFactor = 0.5;  // for RI model
			Par.ZFractionYieldToEject = 0.;
			Par.ZFractionYieldToHot = 0.;
		}
	}

	Par.StarFormationCriticalSurfaceDensity = 10; //Msun/pc^2 no use for the Krumholz model
	Par.StarFormationEfficiency = 0.0; // no use for Krumholz model
	Par.EntropyRatio = 1;

	init_cosmo();

	read_cooling_function();

	init_random_number_gstructure();

	if (Mah_simu) N_halo = read_simu_mah();

	printf("check in Main: Mah_sim=%d Metal_gas_evolu=%d N_halo=%d\n", Mah_simu, Metal_gas_evolu, N_halo);
	init_file();

	return 0;
}

int finalize_run(void)
{
	free_cooling_table();
	if(Mah_simu) free_simu_mah();
	free_random_number_gstructure();
	close_file();
	return 0;
}


int run_galaxy(double *params, int nparams, double *preds, int npreds, int mode, int irun)
{
	int ihalo, ibuf;
	struct galaxy gal;
	char fname_pred[200];

	set_varying_parameters(params, nparams);

	if(Write_pred_file && Write_pred_saparately)
	{
		sprintf(fname_pred, "pred_%04d_z%3.1f.dat", irun, Redshift_end);
		fp_pred=fopen(fname_pred, "w");
	}

	for (ihalo=0; ihalo<N_halo; ihalo++)
	{
		if (Mah_simu) select_simu_mah(ihalo);
		else Mass_bin = LogHaloMassArray[ihalo];

		init(&gal);

		evolve_galaxy(&gal, mode);
		printf("working on %d / %d halo logMass=%g\n", ihalo, N_halo, Mass_bin );

		if(Write_pred_file)
			fprintf(fp_pred, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
			gal.MassHalo, gal.RadiusHalo*1e3, gal.VelocityVirial, gal.MassStar, gal.RadiusHalfStar*1e3, gal.MassCold, gal.MassColdAtomic, gal.RateStarFormation, gal.MassHot, gal.MassEject, gal.MetalHot, gal.MetalCold, gal.MetalStar, gal.MetalEject);

		printf("Mvir=%g M*=%g R*=%g Mc=%g Mh=%g S_vir=%g S_ph=%g f_spin=%g\n", gal.MassHalo, gal.MassStar, gal.RadiusHalfStar*1e3,gal.MassCold, gal.MassHot, gal.EntropyVirial, Par.PreheatEntropy, Par.DiskRadiusFactor);

		set_predictions(preds, npreds, &gal, ihalo);
	}

	if(Write_pred_file && Write_pred_saparately) fclose(fp_pred);
	return 0;
}

void init_file(void)
{
	char fname_pred[200];
	char fname_hist[200];
	char fname_disc[200];

	sprintf(fname_pred, "sample_z%3.1f.dat", Redshift_end);
	fp_list=fopen("list.dat","w");
	fp_pred=fopen(fname_pred, "w");


	{
		sprintf(fname_hist, "hist.dat");
		fp_hist=fopen(fname_hist,"w");
		fprintf(fp_hist, "#z thubble MassHalo MassHot MassCold MassStar MassEject MassColdAtomic MassColdMolecular MassColdIonized RadiusHalo RadiusDisc RadiusHalfCold RadiusHalfStar RadiusCooling ConcenHalo RateHaloAccretion RateCooling RateStarFormation RateOutflow VelocityVirial EntropyVirial TimeCooling MetalHot MetalCold MetalStar MetalEject MassBin\n");
	}

	{
		sprintf(fname_disc, "disc.dat");
		fp_disc=fopen(fname_disc,"w");
		fprintf(fp_disc, "#i RadiusInner SDensityCold SDensityStar SDensityColdMolecular SDensityColdAtomic RadiusOuter MassProfHalo MassProfStar MassProfCold MassProfHot DensityProfHot TemperatureProfHot CoolingRate CoolingTime SDensitySFR SDensityOFR MassProfDM MassProfDMContracted MassMetalCold MassMetalStar SDensityMetalCold SDensityMetalStar MetallicityCold MetallicityStar MassBin RadiusIso\n");
	}

	//printf("done init_file\n");
}

void close_file(void)
{
	fclose(fp_pred);
	fclose(fp_list);
	fclose(fp_hist);
	fclose(fp_disc);
}
