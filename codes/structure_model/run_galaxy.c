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
int Do_preheating = 1;
int Star_formation_model = 1;
int Do_reinfall = 0;
int N_halo = 11;
float Mass_bin = 1;//12.0;
float LogHaloMassArray[11]={10.0,10.25, 10.5, 10.75, 11.0, 11.25, 11.5, 11.75, 12.0,12.25,12.5};
int Resize_radius_bins=1;
int Write_pred_file=1;
int Write_pred_saparately=0;
int Write_hist_file=1;
int Write_SFH_file=1;
int Write_ZFH_file=1;
int Write_prof_file=1;
int Write_snap_file=1;
double Redshift;
double Redshift_end=0.0;

char Mah_file_name[400];
struct parameter Par;

FILE *fp_head;
FILE *fp_hist;
FILE *fp_SFH;
FILE *fp_ZFH;
FILE *fp_pred;
FILE *fp_disc;
FILE *fp_list;
FILE *fp_snap;

int setup_run(char *fname)
{

	UnitMass_in_g = MSUN_IN_G;
	UnitLength_in_cm = PC_IN_CM*1e6;
	UnitTime_in_Megayears = 1e3 ;
	UnitTime_in_Second = 1e9*365*24*3600.;
	UnitVelocity_in_cm_per_s = 1e5;

	read_parameter_file(fname);

	init_parameters();

	init_cosmo();

	read_cooling_function();

	if (Mah_simu) N_halo = read_simu_mah();

	printf("check in Main: Mah_sim=%d Metal_gas_evolu=%d N_halo=%d\n", Mah_simu, Metal_gas_evolu, N_halo);

	return 0;
}

int finalize_run(void)
{
	free_cooling_table();
	if(Mah_simu) free_simu_mah();

	return 0;
}


int run_galaxy(double *params, int nparams, double *preds, int npreds, int mode, int irun)
{
	int ihalo;
	struct galaxy gal;
	char fname_pred[200];

	init_file(irun);
    fprintf(fp_head, "%d %d\n", N_halo, N_RADIUS_BIN);
    fprintf(fp_head, "#ntbin\n");

	set_varying_parameters(params, nparams);

	init_random_number_gstructure();

	if(Write_pred_file && Write_pred_saparately)
	{
		sprintf(fname_pred, "%s/pred_%04d_z%3.1f.dat", OutputDir, irun, Redshift_end);
		fp_pred=fopen(fname_pred, "w");
	}

	for (ihalo=0; ihalo<N_halo; ihalo++)
	{
		if (Mah_simu) select_simu_mah(ihalo);
		//else Mass_bin = LogHaloMassArray[ihalo];

		init(&gal);

		evolve_galaxy(&gal, mode);
		printf("working on %d / %d halo logMass=%g\n", ihalo, N_halo, Mass_bin );

		if(Write_pred_file)
			fprintf(fp_pred, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
			gal.MassHalo, gal.RadiusHalo*1e3, gal.VelocityVirial, gal.MassStar, gal.RadiusHalfStar*1e3, gal.MassCold, gal.MassColdAtomic, gal.RateStarFormation, gal.MassHot, gal.MassEject, gal.MetalHot, gal.MetalCold, gal.MetalStar, gal.MetalEject);

		printf("z=%g Mvir=%g M*=%g R*=%g Mc=%g Mh=%g S_vir=%g S_ph=%g f_spin=%g\n", gal.z, gal.MassHalo, gal.MassStar, gal.RadiusHalfStar*1e3,gal.MassCold, gal.MassHot, gal.EntropyVirial, Par.PreheatEntropy, Par.DiskRadiusFactor);

		set_predictions(preds, npreds, &gal, ihalo);
	}

	if(Write_pred_file && Write_pred_saparately) fclose(fp_pred);

	close_file();
	free_random_number_gstructure();
	return 0;
}

void init_file(int irun)
{
	char fname_head[200];
    char fname_pred[200];
	char fname_hist[200];
    char fname_SFH[200];
    char fname_ZFH[200];
	char fname_disc[200];
    char fname_snap[200];

<<<<<<< HEAD
    {
        sprintf(fname_head, "header.dat");
        fp_head=fopen(fname_head,"w");
        fprintf(fp_head, "#nMAH nrbin\n");
    }
    
	sprintf(fname_pred, "sample_z%3.1f_m%d.dat", Redshift_end, irun);
=======
	sprintf(fname_pred, "%s/sample_z%3.1f_m%d.dat", OutputDir, Redshift_end, irun);
>>>>>>> origin/adding-metallicity-profile
	fp_pred=fopen(fname_pred, "w");


	{
		sprintf(fname_hist, "%s/hist.dat", OutputDir);
		fp_hist=fopen(fname_hist,"w");
		fprintf(fp_hist, "#z thubble MassHalo MassHot MassCold MassStar MassEject MassColdAtomic MassColdMolecular MassColdIonized RadiusHalo RadiusDisc RadiusHalfCold RadiusHalfStar RadiusCooling ConcenHalo RateHaloAccretion RateCooling RateStarFormation RateOutflow VelocityVirial EntropyVirial TimeCooling MetalHot MetalCold MetalStar MetalEject \n");
	}
    
    {
        sprintf(fname_SFH, "SFH.dat");
        fp_SFH=fopen(fname_SFH,"w");
    }
    
    {
        sprintf(fname_ZFH, "ZFH.dat");
        fp_ZFH=fopen(fname_ZFH,"w");
    }

	{
		sprintf(fname_disc, "%s/disc.dat", OutputDir);
		fp_disc=fopen(fname_disc,"w");
		fprintf(fp_disc, "#i RadiusInner RadiusOuter RadiusIso SDensityCold SDensityStar SDensityColdMolecular SDensityColdAtomic MassProfHalo MassProfStar MassProfCold MassProfHot DensityProfHot TemperatureProfHot CoolingRate CoolingTime SDensitySFR SDensityOFR SDensityCAR StellarAge MassProfDM MassProfDMContracted MassMetalCold MassMetalStar SDensityMetalCold SDensityMetalStar MetallicityCold MetallicityStar MassStar MassHalo\n");
	}

    {
        sprintf(fname_snap, "%s/snap.dat", OutputDir);
        fp_snap=fopen(fname_snap,"w");
        fprintf(fp_snap, "#z i RadiusInner RadiusOuter RadiusIso SDensityCold SDensityStar SDensityColdMolecular SDensityColdAtomic SDensitySFR  SDensityOFR SDensityCAR StellarAge MetallicityCold MetallicityStar MassProfStar MassStar MassHalo\n");
    }
 
    

	//printf("done init_file\n");
}

void close_file(void)
{
	fclose(fp_head);
    fclose(fp_pred);
	fclose(fp_hist);
    fclose(fp_SFH);
    fclose(fp_ZFH);
	fclose(fp_disc);
    fclose(fp_snap);
}
