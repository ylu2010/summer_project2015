#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "variables.h"
#include "proto.h"
#include "cosmo.h"

#define PROTONMASS 1.6726e-24
#define BOLTZMANN  1.3806e-16


static float metallicity[]={
    -5.00,
    -3.00,
    -2.00,
    -1.50,
    -1.00,
    -0.05,
    -0.00,
    +0.05};
static float LogLambda[8][91];
static float *ccore, *mratio;

void read_cooling_function(void)
{
    FILE *fd;
    int i,j;
    float nt, ne, mu;
    float logT, temp, loglambda_norm, loglambda_net;
    char filename[256];
    char *coolfilename[8]={
            "mzero.cie",
            "m-30.cie",
            "m-20.cie",
            "m-15.cie",
            "m-10.cie",
            "m-05.cie",
            "m-00.cie",
            "m+05.cie"          };

//      LogLambda = matrix(0,7,1,91);
    for(i=0; i<8; i++)metallicity[i]+=log10(SolarMetallicity);
    for(i=0; i<8; i++)
    {
        sprintf(filename,"/Users/Nathan/project_data/ModelTables_cb07/CoolFunctions/stripped_%s", coolfilename[i]);
        if(!(  fd=fopen(filename,"r")))
        {
            printf("I can not open the file '%s'.",filename);
            exit(1);
        }

        for(j=0; j<91; j++)
        {
            fscanf(fd, "%f %f %f %f %f %f %f %f %f %f %f %f",
                &logT, &ne, &temp, &nt, &loglambda_net, &loglambda_norm,&temp, &temp, &temp, &temp, &temp, &mu);
            LogLambda[i][j]=loglambda_norm +log10(nt*ne/(nt+ne)/(nt+ne)/mu);
        }
        fclose(fd);
    }
}



float get_LogLambda(float logtemp, float logmetal)
{
	int Z,j;
  	float Loglambda;

  	if(!(Metal_gas_evolu))logmetal=log10(SolarMetallicity*MinimumMetallicityRelativeToSolar);
  	if(logmetal<metallicity[0])logmetal=metallicity[0];
  	if(logmetal>metallicity[7])logmetal=metallicity[7];
  	Z=0;
  	while(logmetal>metallicity[Z+1])
    {
      		Z++;
    }
  	if(logtemp<4.00)logtemp=4.00;
  	else if(logtemp>8.50)logtemp=8.50;
  	j=(int)((logtemp-4.00)/0.05)+1;
  	Loglambda = LogLambda[Z][j]+(LogLambda[Z+1][j]-LogLambda[Z][j])/(metallicity[Z+1]-metallicity[Z])*(logmetal-metallicity[Z]);
  	Loglambda+= (LogLambda[Z][j+1]-LogLambda[Z][j])/0.05*(logtemp-(4.00+(j-1)*0.05));
  	return Loglambda;
}



static float rcool_old = 0.0;

double cooling_rate(struct galaxy *gal, double hubble_time, double dt, int model)
{
	double rcore, arv, arc, rho_0;
	double t_gas, z_gas, gas;
	double logLambda, cof;
	double cool_time, rcool, rcool1, rcool2, dhotgas, cool_rate, dyn_time;
	float mr;
	float con, dcon;
	int i;
	unsigned long j;

	double H0kpc=0.1, H0mpc=100.;
	double mvir, rvir, vvir, hotgas, coldgas;

	mvir = gal->MassHalo;
	rvir = gal->RadiusHalo;
	hotgas = gal->MassHot;
	coldgas = gal->MassCold;
	vvir = gal->VelocityVirial;

    t_gas = 35.9 * vvir*vvir; //vvir in km/s
	z_gas = 1e-1;
	dyn_time = rvir*UnitLength_in_cm*1e-5/vvir/UnitTime_in_Second;
	//printf("cooling: %g %g %g\n", mvir, t_gas, dyn_time);
   	if(hotgas<=0.0)
	{
  		cool_rate = 0.;
	}
	else
	{
		logLambda=get_LogLambda(log10(t_gas), log10(z_gas));
		//cof=1./0.28086 *PROTONMASS*BOLTZMANN*t_gas/pow(10,logLambda)
		//	/UnitMass_in_g*(UnitLength_in_cm*UnitLength_in_cm*UnitLength_in_cm);
		cof= 0.9 *PROTONMASS*BOLTZMANN*t_gas/pow(10,logLambda)
			/UnitMass_in_g*(UnitLength_in_cm*UnitLength_in_cm*UnitLength_in_cm);
		cof*=1./UnitTime_in_Second;//H0mpc*UnitVelocity_in_cm_per_s/UnitLength_in_cm;
		switch(model)
        {
            case 1: // croton
                rcore = 0.0;
                rho_0 = hotgas/4./PI/rvir;
                cool_time = dyn_time;
                rcool = rho_0*cool_time/cof;
                rcool = sqrt(rcool);
                cool_rate = hotgas*sqrt(rcool*rcool-rcore*rcore)/(2.*rvir*cool_time);
                //printf("cooling1: %g %g %g %g\n", rcool/rvir, rvir, hotgas, cool_rate);
                if(rcool > rvir) 
				{
					cool_rate = 0.5*hotgas/dyn_time;
                	//cool_rate = hotgas/dt;
					rcool = rvir;
				}
                gal->RadiusCooling = rcool;
                //printf("cooling2: %g %g %g %g\n", gal->RadiusCooling, rvir, hotgas, cool_rate);
                break;
            case 2:  // kang
                rcore = 0.0;
                rho_0 = hotgas/4./PI/rvir;
                cool_time = hubble_time;
                rcool = rho_0*cool_time/cof;
                rcool = sqrt(rcool);
                cool_rate = hotgas*sqrt(rcool*rcool-rcore*rcore)/(2.*rvir*cool_time);

                if(rcool > rvir) 
				{
					cool_rate = hotgas/cool_time;
					rcool = rvir;

				}
                gal->RadiusCooling = rcool;
                printf("cooling2: %g %g %g %g\n", gal->RadiusCooling, rvir, hotgas, cool_rate);
                break;
            default:
            	printf("You did not choose an appropriate model!\n");
            	printf("exit from cooling!\n");
            	exit(0);
        }

	   	if(cool_rate*dt>hotgas) cool_rate=hotgas/dt;
	}
	return cool_rate;
}


void free_cooling_table(void)
{
//	free_matrix(LogLambda, 0,7,1,91);
}

double specific_angular_momentum_profile(double r)
{
	double j;

	j = pow(r, 1.1);
	return j;
}

double cooling_rate_shell(struct galaxy *gal, double hubble_time, double dt, int model)
{
	double rcore, arv, arc, rho_0;
	double t_gas, z_gas, logz_gas, gas;
	double logLambda, cof, fac, den;
	double cool_time, rcool, rcool1, rcool2, dhotgas, cool_rate, dyn_time, ff_time, cr=0;
	double j_specific=0, j_gas=0, j_halo=0, halo_mass, mstar, hr;
	float mr;
	float con, dcon;
	int i;
	unsigned long j;

	double H0kpc=0.1, H0mpc=100.;
	double mvir, rvir, vvir, hotgas, coldgas, rdisk;

	mstar = gal->MassStar;
	mvir = gal->MassHalo;
	hotgas = gal->MassHot;
	rvir = gal->RadiusHalo;
	vvir = gal->VelocityVirial;
	rdisk = gal->RadiusDisc;
	//hr = 7e-5*mstar * vvir*vvir;
	//hr = 7e-5*mstar * vvir*130;
	//hr = 1.5* mstar;
	//hr = 1.0e-10*mstar *mstar;
	//hr = 0.0;
	//hr = mstar * mstar/rdisk * 0.005*1e-11;
	//hr = 7e-5*mstar/rdisk/rdisk * 130 * 130 * 0.005 * 0.005;
	//printf("t=%g rdisk:%g\n", hubble_time, rdisk);

    t_gas = 35.9 * vvir*vvir; //vvir in km/s
	//logz_gas = log10(0.02*1.e-1);
    logz_gas = log10(gal->MetalHot/gal->MassHot);

	/* tests for Peter Mitchell start:
	logz_gas = log10(0.02);

	for(i=0; i<gal->nbin; i++)
	{
		gal->TemperatureProfHot[i] = 6e5;
		gal->DensityProfHot[i] = 1e3*1e9;
	}

	// tests for Peter Mitchell end */

	dyn_time = rvir*UnitLength_in_cm*1e-5/vvir/UnitTime_in_Second;
	fac = 1.5e-24*BOLTZMANN
			/UnitMass_in_g*(UnitLength_in_cm*UnitLength_in_cm*UnitLength_in_cm) / UnitTime_in_Second;
	//printf("cooling: %g %g %g %g\n", mvir, t_gas, dyn_time, fac);
	rcool=0;
	gal->SpinCooling = 1.0;
   	if(hotgas<=0.0)
	{
  		cool_rate = 0.;
  		rcool = 0.0;
  		cool_time = 99.99;
  		cr = 0.0;
	}
	else
	{
		j_gas = 0.0; j_halo=0.0; halo_mass = 0.0;
		for(i=0, hotgas = 0, cr=0; i<gal->nbin; i++)
		{
			if(gal->RadiusOuter[i] >= gal->RadiusHalo) break;

			ff_time = gal->RadiusOuter[i]/rvir * dyn_time;
			t_gas = gal->TemperatureProfHot[i];
			den = gal->DensityProfHot[i];
			logLambda=get_LogLambda(log10(t_gas), logz_gas);
		//cof=1./0.28086 *PROTONMASS*BOLTZMANN*t_gas/pow(10,logLambda)
		//	/UnitMass_in_g*(UnitLength_in_cm*UnitLength_in_cm*UnitLength_in_cm);
			cof = t_gas/den/pow(10,logLambda);
			cool_time = cof * fac; // cooling timescale
//printf("test for Peter Mitchel: T=%g Z-Zsun=%g rho=%g tcool=%g t_ff=%g\n", t_gas, logz_gas-log10(0.02), den, cool_time, ff_time);
			cool_rate = (gal->MassProfHot[i] - hotgas) / dmax(cool_time, ff_time);
			//cool_rate = (gal->MassProfHot[i] - hotgas) / cool_time;
			//cool_rate = (gal->MassProfHalo[i]-halo_mass)*1e-6;
			cr += cool_rate;
			gal->CoolingRate[i] = cool_rate;
			j_specific = specific_angular_momentum_profile(gal->RadiusOuter[i]);
			j_gas +=  j_specific * cool_rate;
			j_halo +=  j_specific * (gal->MassProfHalo[i]-halo_mass);
                //printf("cooling1: %d %g %g %g %g %g %g\n", i, gal->RadiusOuter[i], gal->DensityProfHot[i], gal->TemperatureProfHot[i], cool_time, ff_time, dt);
            if(cool_time > dyn_time)
			{
				rcool = gal->RadiusOuter[i];
			}

            //printf("cooling2: %g %g %g %g\n", gal->RadiusOuter[i], rvir, hotgas, cool_rate);
            hotgas = gal->MassProfHot[i];
            halo_mass = gal->MassProfHalo[i];
            if (i==0) gal->TimeCooling = cool_time;
            gal->CoolingTime[i] = cool_time;
		}
		if(cr > 0.0) gal->SpinCooling = j_gas/j_halo*halo_mass/cr;

	}
   	///exit(0); //test for Peter Mitchel
   	gal->TimeCooling = cool_time;
   	//gal->TimeCooling = t_gas;
   	//if(cr > 0.0) gal->SpinCooling = j_gas/j_halo*halo_mass/cr;
   	//printf("cooling: %g %g %g %g %g %g %g\n", gal->MassHalo, t_gas, den*UnitMass_in_g/(UnitLength_in_cm*UnitLength_in_cm*UnitLength_in_cm)/PROTONMASS, cool_time, gal->TimeCooling, hubble_time, UnitTime_in_Second/1e9/3.15360e+07);
   	//printf("cooling: %g %g j=%g %g %g %g %g %g\n", hubble_time, gal->MassHot, j_gas, j_halo, j_gas/j_halo, cr/halo_mass, gal->SpinCooling, gal->TimeCooling);
   	//printf("spin in cooling: %g %g\n", gal->z, j_gas/j_halo*halo_mass/cr);
   	gal->RadiusCooling = rcool;
//printf("cooling: z=%g hr=%g cr=%g\n", gal->z, hr, cr);
	//if (Do_preheating) cr = dmax(cr - hr, 0.0);
   	return cr;
}

