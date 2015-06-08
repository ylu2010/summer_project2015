/*
 * evolve_galaxy.c
 *
 *  Created on: Jun 14, 2012
 *      Author: luyu
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "variables.h"
#include "proto.h"
#include "cosmo.h"

void reset_parameter(double z, double m)
{
	//if ( m/mah(0.0) > 1./10)
	//printf("c=%g\n", halo_concentration_prada(0, mah(0.0)));
	if (m/mah(0.0) > 0.01/halo_concentration_prada(0, mah(0.0)))
	{
		Do_preheating = 1;
		Par.PreheatEntropySlope = 0.2;
		Par.PreheatEntropy = 12;//17;//12;//* pow(pow(10.,Mass_bin)/1e12, Par.PreheatEntropySlope);
		Par.EntropyProfileIndex = 0.0;
		Par.DiskRadiusFactor = 0.7;
		Par.ZFractionYieldToEject = 0.;// for PR model
		Par.ZFractionYieldToHot = 0.2;  // for PR model
		Par.GalaxyHeatingEfficiency = 1.0;//1.;//1.1;
		Par.SNLoadingFactor= 1.;
		Par.SNLoadingFactorIndex = 0.;
	}
	else
	{
		Do_preheating = 0;
		Par.PreheatEntropySlope = 0.;
		Par.PreheatEntropy = 0.0;
		Par.EntropyProfileIndex = 1.1;//4./3;
		Par.DiskRadiusFactor = 0.8;    // for EJ model
		Par.ZFractionYieldToEject = 0.;// for EJ model
		Par.ZFractionYieldToHot = 0.2;  // for EJ model
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
}

void evolve_galaxy(struct galaxy *gal, int mode)
{
	double z_min, z_max, mvir_in;
	double z0, dz, dz0, z;
	double mh, ms, mc, mhot;
	double dmhdt, thubble, dt, dmh;
	double rho, v, dar, bar, f_hot_accretion, f_cloud_accretion, factor;


	factor = 1.e3 * ELECTRONVOLT * pow(BaryonFrac/0.6/PROTONMASS*UnitMass_in_g, 2./3)/BOLTZMANN / (UnitLength_in_cm * UnitLength_in_cm);

	z_min = z0 = Redshift_end; z_max = 10;
	mvir_in = mah(z_max);
 	dz = (z_max-z_min)/400;
	dz0 = dz;

	Redshift = z = z_max;
	gal->z = z;
	gal->MassHalo = mh = mvir_in;
	gal->MassStar = ms = 0;
	gal->MassCold = mc = 0;
	dz = dz0;

	//reset_parameter(z, mh); // *******************///

	rho = Delta_vir(z) * rho_crit(z) * xhubble * xhubble;
	gal->VelocityVirial = v = sqrt(G) * pow(4./3*M_PI*rho, 1./6) * pow(mh, 1./3);
	//gal->VelcoityVirial = v = pow(G*mh*xH(z) * sqrt(0.5*Delta_vir(z)), 1./3);  // just another way

	gal->RadiusHalo = pow(3.*mh/(4.*M_PI*rho), 1./3);
	gal->ConcenHalo = halo_concentration(z, gal);
	gal->VelocityMax = halo_vmax(gal->VelocityVirial, gal->ConcenHalo);
	gal->TemperatureVirial = 35.9 * v * v;
	gal->EntropyVirial = gal->TemperatureVirial/pow(rho, 2./3) / factor;
	gal->MassHot = mhot = mh * hot_gas_accretion_fraction(gal, z);
	gal->RadiusHalfStar = 0.015 * gal->RadiusHalo;
	gal->RadiusHalfCold = gal->RadiusHalfStar * 2.6;
	gal->RadiusDisc = disk_radius(gal);
	disc_mass_composition(gal);

	while (z >= z_min-0.5*dz)
	{
		thubble = hubble_time(z) * xH0_recip;
		dt = thubble - hubble_time(z+dz) * xH0_recip;
		dmhdt = halo_mass_accretion_rate( mh, z);
		gal->RateHaloAccretion = dmhdt;
		dmh = dmhdt * dt;

		//reset_parameter(z, mh);

		if(Mah_simu == 0)
		{
			if( dmh > 0.01*mh ) {	// if too fast, reduce timestep
				dz *= 0.5;
				continue;
			}
			if( dmh < 0.001*mh ) {	// if too slow, increase dz
				dz *= 2;
				continue;
			}
			//dz = dz0;
		}
		//else dz = dz0;
		else dz = z_hubble_time(thubble) - z_hubble_time(thubble + 13.6/400);
		//printf("z=%g dz=%g\n",z, dz);
		rho = Delta_vir(z) * rho_crit(z) * xhubble * xhubble;
		gal->VelocityVirial = v = sqrt(G) * pow(4./3*M_PI*rho, 1./6) * pow(mh, 1./3);
		//gal->VelcoityVirial = v = pow(G*mh*xH(z) * sqrt(0.5*Delta_vir(z)), 1./3);  // just another way
		gal->RadiusHalo = pow(3.*mh/(4.*M_PI*rho), 1./3);
		gal->ConcenHalo = halo_concentration(z, gal);
		gal->TemperatureVirial = 35.9 * v * v;
		gal->EntropyVirial = gal->TemperatureVirial/pow(rho, 2./3) / factor;
		gal->RadiusDisc = disk_radius(gal); //Par.DiskRadiusFactor *(1.+z)/2 * 0.035/sqrt(2.) * gal->RadiusHalo;

		halo_mass_profile(gal);

		//hot_gas_profile(gal);
		//hot_gas_profile_isenthermal(gal);
		//gal->ConcenHalo = 10.;
		//gal->MassHot = BaryonFrac * gal->MassHalo;
		hot_gas_profile_power_law_entropy(gal);

		//cold_gas_accretion(gal, thubble, dt);
		cold_gas_accretion_surface(gal, thubble, dt);

		disc_mass_composition(gal);
		//star_formation_global(gal, dt);
		//star_formation_surface(gal, thubble+0.5*dt, dt);

		if(Do_reinfall)
		{
			star_formation_surface_molecule_with_guo2011_feedback(gal, thubble+0.5*dt, dt); // a model mimicking the feedback model of Guo et al. 2011
			ejected_gas_reincorporation(gal, thubble+0.5*dt, dt);
		}
		else
			star_formation_surface_molecule(gal, thubble+0.5*dt, dt); // the model used in the preheating paper

		////star_formation_surface_limit(gal, dt);

		//disc_mass_composition(gal);
		//printf("%g %g %g %g %g %g %g %g\n", z, thubble, gal->MassHalo, gal->MassHot, gal->MassCold, gal->MassStar, gal->RadiusHalo, gal->VelocityVirial);

		//if(write_in_file) fprintf(fp_halo, "%g %g %g %g %g %g %g %g %g %g %g %g %g\n", z, thubble, gal->MassHalo, gal->MassHot, gal->MassCold, gal->MassStar, gal->RateHaloAccretion, gal->RateStarFormation,
		//		gal->MassEject, gal->MassColdMolecular, gal->MassColdAtomic, gal->MassColdIonized, gal->RateCooling);


		//if(gal->MassStar > 0.00001 * gal->MassHalo && gal->MassCold > 0.00001 * gal->MassHalo)
		{
			gal->RadiusHalfStar = interpolate_bipoint(gal->MassProfStar, gal->RadiusOuter, N_RADIUS_BIN, 0.5*gal->MassProfStar[N_RADIUS_BIN-1], 0);  // in Mpc
			gal->RadiusHalfCold = interpolate_bipoint(gal->MassProfCold, gal->RadiusOuter, N_RADIUS_BIN, 0.5*gal->MassProfCold[N_RADIUS_BIN-1], 0);  // in Mpc
		}
		/*
		else
		{
			gal->RadiusHalfStar = 0.015 * gal->RadiusHalo;
			gal->RadiusHalfCold = gal->RadiusHalfStar * 2.6;
		}
		*/
		if(Write_hist_file)
		{
			f_hot_accretion = hot_gas_accretion_fraction(gal, z);

			gal->MassHot = dmax(gal->MassHalo * f_hot_accretion - (gal->MassStar+gal->MassCold+gal->MassEject), 0);
			disc_mass_composition(gal);
			//printf("mcold=%g matom=%g %g\n", gal->MassCold, gal->MassColdAtomic, f_hot_accretion);
			fprintf(fp_hist, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n", gal->z, thubble,
				gal->MassHalo, gal->MassHot, gal->MassCold, gal->MassStar, gal->MassEject,
				gal->MassColdAtomic, gal->MassColdMolecular, gal->MassColdIonized,
				gal->RadiusHalo, gal->RadiusDisc, gal->RadiusHalfCold, gal->RadiusHalfStar, gal->RadiusCooling, gal->ConcenHalo,
				gal->RateHaloAccretion, gal->RateCooling, gal->RateStarFormation, gal->RateOutflow,
				gal->VelocityVirial, gal->EntropyVirial, gal->TimeCooling, gal->MetalHot, gal->MetalCold, gal->MetalStar, gal->MetalEject);
		}
		z -= dz;
		Redshift = z;
		mh += dmh;
		dmhdt = dar;
		gal->z = z;
		gal->MassHalo = mh;

		f_hot_accretion = hot_gas_accretion_fraction(gal, z);

		gal->MassHot = dmax(gal->MassHalo * f_hot_accretion - (gal->MassStar+gal->MassCold+gal->MassEject), 0);
		f_cloud_accretion = halo_cold_gas_accretion_fraction(gal, z);
		gal->MassCloud += f_cloud_accretion * dmh;
		//printf("evolve: z=%g %g %g %g %g %g\n", z, mh, gal->VelocityVirial, v, f_hot_accretion, gal->MassStar);
		//gal->MassHot = mhot;

		//Splitting Metals evenly into radius bins
		int i;
		for(i=0; i<gal->nbin; i++)
		  {
		    gal->MassMetalCold[i] = (gal->MetalCold)/(gal->nbin);
		    gal->MassMetalStar[i] = (gal->MetalStar)/(gal->nbin);
		    gal->SDensityMetalCold[i] = (gal->MassMetalCold[i])/(2. * M_PI * (((gal->RadiusOuter[i])*(gal->RadiusOuter[i])) - ((gal->RadiusInner[i])*(gal->RadiusInner[i]))));
		    gal->SDensityMetalStar[i] = (gal->MassMetalStar[i])/(2. * M_PI * (((gal->RadiusOuter[i])*(gal->RadiusOuter[i])) - ((gal->RadiusInner[i])*(gal->RadiusInner[i]))));
		  }

	}
	//printf("evolve:%g\n",gal->RateStarFormation );
	//adiabatic_contraction(gal);
	if(Write_prof_file) print_galaxy(gal);
}

void print_galaxy(struct galaxy *gal)
{
	int i;
	for(i=0; i<gal->nbin; i++)
	{
		fprintf(fp_disc, "%d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
				i, (gal->RadiusInner[i])*1e3, (gal->SDensityCold[i])/1e12, (gal->SDensityStar[i])/1e12, (gal->SDensityColdMolecular[i]/1e12), (gal->SDensityColdAtomic[i]/1e12),
				(gal->RadiusOuter[i])*1e3, gal->MassProfHalo[i], gal->MassProfStar[i], gal->MassProfCold[i], gal->MassProfHot[i],
				gal->DensityProfHot[i]/1e9, gal->TemperatureProfHot[i], gal->CoolingRate[i],gal->CoolingTime[i],gal->SDensitySFR[i],
			        gal->MassProfDM[i], gal->MassProfDMContracted[i], gal->MassMetalCold[i], gal->MassMetalStar[i], gal->SDensityMetalCold[i]/1e12, gal->SDensityMetalStar[i]/1e12);
	}
}
