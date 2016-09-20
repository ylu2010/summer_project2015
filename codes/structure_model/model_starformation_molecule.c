/*
 * model_starformation_molecule.c
 *
 *  Created on: Jun 14, 2012
 *      Author: luyu
 *      To do: add stellar mass return
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>

#include "variables.h"
#include "proto.h"
#include "cosmo.h"

double stellar_mass_loss(double t)
{
	double c0, t0, tau, fml=0;
	/* Leitner paper
	c0 = 0.046;
	tau = 2.76e-4; // Gyears

	fml = c0 * log( t / tau +1);
	*/
	c0 = 0.05;
	t0 = 3e-3; // Gyr
	tau = 3.76e-4; // Gyr
	if(t > t0) fml = c0 * log( (t-t0) / tau +1);
	//printf("fml: %g %g\n", t, fml);

	return fml;
}

double molecular_fraction_krumholz(double sig_gcm, double z0)
{
    double clumping = 5.0;
    double den, tauc, x, s, f_h2=0.0;

    clumping = 5+ 5*(0.5*(1+gsl_sf_erf((Redshift - 1.2)/0.5)));
    tauc = 0.066 * clumping *z0 * sig_gcm * 1e-12;
    x = 3.1 * (1.+3.1*pow(z0, 0.365))/4.1;
    s = log(1.+ 0.6*x + 0.01*x*x )/(0.6*tauc);

    if (s < 2 ) f_h2 = 1 - 0.75 * (s/(1.+0.25*s));
    //f_h2 = 1 - 0.75 * (s/(1.+0.25*s));
//printf("mf_krumhoz: %g %g\n", sig_gcm, f_h2);
    return f_h2;
}

double molecular_fraction_fu_blitz(double sig_cg, double sig_star, double sig_star0)
// need test !!!! in unit of Msun/pc^2
{
	double f;

	f = 1.38e-3 * pow(sig_cg*sig_cg + 0.1 * sig_cg * sqrt(sig_star * sig_star0), 0.92);

	return dmin(f, 1);
}

double molecular_fraction(double den, double z0)
{
	double f=0;

	// simple cut-off model
	//if (den >= Par.StarFormationCriticalSurfaceDensity*1e12) f = 1;
	//else f = 0;

	// Krumholz model
	f = molecular_fraction_krumholz(den, z0);
	return f;
}


double star_formation_efficiency_krumholz(double sig_gcm)
{
    double sig_gcm0 = 85.e12, t0, fff, e0; // 85Msun/pc^2

    t0 = pow(sig_gcm/sig_gcm0, 0.33);
    e0 = 1./2.6;
    //e0 = Par.StarFormationEfficiency;
    if(sig_gcm >= sig_gcm0) fff = t0* e0;
    else fff = 1./t0 * e0;
    return fff;
}

double star_formation_efficiency_kennicutt(double den)
{
	double a = 0.167; // for the efficiency for gas density in units of Msun/Gyr/Mpc^2
	double n = 0.4;
	double sigma = den * 1e-12;  // convert from Msun/Mpc^2 to Msun/pc^2
	double eff = a * pow(sigma, n);

	return eff;
}
double star_formation_efficiency(double den, double tau)
{
	double sfe;
	// simple model
	//sfe = Par.StarFormationEfficiency / tau * pow((den/(Par.StarFormationCriticalSurfaceDensity*1e12)), 0.4);

	// Kennicutt model
	//sfe = star_formation_efficiency_kennicutt(den);

	// Krumholz model
	sfe = star_formation_efficiency_krumholz(den);
	//printf("sfe=%g\n", sfe);
	return sfe;
}


void star_formation_global(struct galaxy *gal, double dt)
{
    int i, nbin;
    double z0, u0, f, sig0_cold, sig0_star;
    double reff, r_crit, f_sf, mass_sf, sfr, tdyn;

    nbin = gal->nbin;

    u0 = UnitMass_in_g / (UnitLength_in_cm * UnitLength_in_cm);

    reff = gal->RadiusDisc;
    //tdyn = 1./((gal->VelocityVirial/(gal->RadiusHalo*UnitLength_in_cm/1e5)) * UnitTime_in_Second);
    tdyn=1.;
    r_crit = log(gal->MassCold/(2.*M_PI*reff*reff)) - log(Par.StarFormationCriticalSurfaceDensity * 1e12);
    r_crit = dmax(0., r_crit);
    f_sf = 1.-(1.+ r_crit)*exp(-r_crit);

    if(f_sf > 0.) mass_sf = gal->MassCold * f_sf;
    else mass_sf = 0.0;

    sfr = Par.StarFormationEfficiency * mass_sf / tdyn;

    gal->MassCold -= (1.-Par.StellarMassLossFraction) * sfr * dt;
    gal->MassStar += (1.-Par.StellarMassLossFraction) * sfr * dt;
    gal->MassCold -= Par.SNLoadingFactor * sfr * dt;
    gal->MassEject += Par.SNLoadingFactor * sfr *dt;
    gal->MassColdMolecular = mass_sf - (1.-Par.StellarMassLossFraction+Par.SNLoadingFactor) * sfr *dt;
    gal->MassColdAtomic = gal->MassCold - gal->MassColdMolecular;

    sig0_cold = gal->MassCold/(2.*M_PI * gal->RadiusDisc * gal->RadiusDisc);
    sig0_star = gal->MassStar/(2.*M_PI * gal->RadiusDisc * gal->RadiusDisc);
    
    for(i=0; i<nbin; i++)
    {
    	f = exp(-0.5*(gal->RadiusInner[i]+gal->RadiusOuter[i])/gal->RadiusDisc);
    	gal->SDensityCold[i] = sig0_cold * f;
    	gal->SDensityStar[i] = sig0_star * f;
    }

    gal->RateStarFormation = sfr;
    //printf("%g %g\n", Par.StarFormationEfficiency, Par.SNLoadingFactor );
    //printf("sfr=%g %g %g %g %g %g %g\n", sfr, f_sf, r_crit, reff, gal->RadiusHalo, tdyn, dt);
}

void star_formation_surface(struct galaxy *gal, double t, double dt)
{
    int i, j, k, nbin;
    double z0, u0;
    double dden, area, sfe, sfr, sd_sfr, ofr, sd_ofr, sd_cold_tmp, mstar, mcold, fml;
    double loadingfactor, tdyn, sfe_max=1e33; //0.4e12;

    nbin = gal->nbin;

    u0 = UnitMass_in_g / (UnitLength_in_cm * UnitLength_in_cm);

    //tdyn = 1./((gal->VelocityVirial/(gal->RadiusHalo*UnitLength_in_cm/1e5)) * UnitTime_in_Second);
    tdyn = 1;
    loadingfactor = outflow_massloading_factor(gal);
    //if(gal->z >= 2.01) loadingfactor = 100.;
    for (i=0, sfr=0.0, ofr=0.0, mstar=0.0, mcold=0.0; i<nbin; i++)
    {
    	dden = gal->SDensityCold[i];
    	area = M_PI * (gal->RadiusOuter[i]*gal->RadiusOuter[i] - gal->RadiusInner[i] * gal->RadiusInner[i]);
    	if(dden > 0)
    	{
    		sd_cold_tmp = gal->SDensityCold[i];

    		// model with readjusted star formation efficiency
    		//sfe_max = dden/(gal->SDensityCold[i] * (1+Par.SNLoadingFactor)*dt);
    		// model with readjusted outflow loading factor
    		sfe_max = 1e33;

    		sfe = star_formation_efficiency(sd_cold_tmp, tdyn); // depends on total cold gas density
    		sfe = dmin(sfe_max, sfe);
    		sd_sfr = sfe * dden;

    		j = floor(t/Bin_size_time);
    		gal->SDensitySFH[i][j] += sd_sfr * dt;
    		for(k=0, fml=0; k<j; k++)
    		{
    			fml += gal->SDensitySFH[i][k] *
    					(stellar_mass_loss(t - gal->TimeArray[k]+0.5*dt) -
    							stellar_mass_loss(t - gal->TimeArray[k]-0.5*dt));
    		}
    		//printf("%g %g\n", fml, 0.3*sd_sfr*dt);
    		gal->SDensityStar[i] += (1.-Par.StellarMassLossFraction) * sd_sfr *dt;
    		gal->SDensityStar[i] -= fml;


    		// model with readjusted star formation efficiency
    		/*
    		sd_cold_tmp -= (1.-Par.StellarMassLossFraction) * sd_sfr *dt;
    		sd_ofr = Par.SNLoadingFactor * sd_sfr * dt;
    		sd_cold_tmp -= (1.-Par.StellarMassLossFraction) * sd_ofr;
			*/
    		// model with readjusted outflow loading factor
    		//loadingfactor=dmin(Par.SNLoadingFactor, sd_cold_tmp/(sd_sfr*dt) + Par.StellarMassLossFraction-1.);
    		//sd_ofr = dmin(loadingfactor * sd_sfr * dt, sd_cold_tmp);
    		sd_cold_tmp -= (1.-Par.StellarMassLossFraction) * sd_sfr * dt;
    		sd_cold_tmp += fml;
    		//printf("sf_surface:%d %g %g\n", i, fml, sd_sfr*dt*0.43);

    		sd_ofr = dmin(loadingfactor * sd_sfr * dt, sd_cold_tmp);
    		sd_cold_tmp -= sd_ofr;

    		gal->SDensityCold[i] = sd_cold_tmp;

    		sfr += sd_sfr * area;
    		ofr += sd_ofr * area /dt;
    		/*
    		mstar += gal->SDensityStar[i] * area;
    		mcold += gal->SDensityCold[i] * area;
    		*/
    	}
    	mstar += gal->SDensityStar[i] * area;
    	mcold += gal->SDensityCold[i] * area;
    }
    gal->MassCold = mcold;
    gal->MassStar = mstar;
    gal->MassEject += ofr * dt;

    gal->RateStarFormation = sfr;
}


void disc_mass_composition(struct galaxy *gal)
{
	int i, nbin;
	double area, den, matom, mmole, mion, mstar, zcold, zmetal, fmole;
	double si_ion = 0.0;

//printf("disc_mass_composition: zmetal=%g\n", zmetal);
	nbin = gal->nbin;
	for(i=0, matom=0, mmole=0, mion=0, mstar=0; i<nbin; i++)
	{
		area = M_PI * (gal->RadiusOuter[i]*gal->RadiusOuter[i] - gal->RadiusInner[i] * gal->RadiusInner[i]);
		den = gal->SDensityCold[i] - si_ion * 1e12;
		if (den > 0.0)
		{
		        if ((gal->SDensityCold[i]*area) > 0.0) zcold = gal->MassMetalCold[i]/(gal->SDensityCold[i]*area);
			else zcold = 0.0;

			if (Metal_gas_evolu) zmetal = dmax(zcold / SolarMetallicity, MinimumMetallicityRelativeToSolar);
			else zmetal = 0.5;

			fmole = molecular_fraction( den, zmetal );
/*
			gal->SDensityColdMolecular[i] = gal->SDensityCold[i] * fmole;
			gal->SDensityColdAtomic[i] = gal->SDensityCold[i] * (1-fmole);
*/
			matom += gal->SDensityColdAtomic[i] * area;
			mmole += gal->SDensityColdMolecular[i] * area;
			mion += si_ion * area;
		       
		}
		else
		{
			fmole = 0.0;
			mion += gal->SDensityCold[i] * area;
		}
		gal->SDensityColdMolecular[i] = gal->SDensityCold[i] * fmole;
        gal->SDensityColdAtomic[i] = gal->SDensityCold[i] * (1-fmole);

		mstar += gal->SDensityStar[i] * area;

		gal->MassProfStar[i] = mstar;
		gal->MassProfCold[i] = matom + mmole + mion;
		gal->MassProfNeutral[i] = matom + mmole;
	}
	gal->MassColdAtomic = matom;
	gal->MassColdMolecular = mmole;
	gal->MassColdIonized = mion;
	gal->MassCold = matom + mmole + mion;
}

double outflow_massloading_factor(struct galaxy *gal)
{
	double sigma0, f=0, z, r, p0;
	z= gal->z;
/*
	r = dmin(gal->RadiusHalfStar, gal->RadiusHalo);
	sigma0 = gal->MassStar / (2*M_PI*r*r/1.678/1.678)/1e12;
	if(Do_preheating) f= 1.0*pow(sigma0/800, 6)+Par.SNLoadingFactor * ( 0.5 *  (1.+gsl_sf_erf((z - 2.5)/1.0)));
*/
/*
	r = gal->RadiusOuter[79];
	sigma0 = gal->MassProfStar[79] / (2*M_PI*r*r/1.678/1.678)/1e12;
	*/
	p0 = 0;// * pow(200/gal->VelocityVirial, -1);//Par.SNLoadingFactorIndex);
	//if(Do_preheating) f = (p0+1.0*pow(sigma0/600, 6)) + 1.*(Par.SNLoadingFactor-p0) * ( 0.5 *  (1.+gsl_sf_erf((z - 2.5)/1.0)));
	//if(Do_preheating) f = (p0+0.0*pow(gal->MassStar/2.0e10, 2)) + 1.*(Par.SNLoadingFactor-p0) * ( 0.5 *  (1.+gsl_sf_erf((z - 2.5)/1.0))); //disk paper used this
	if(Do_preheating) //f = (p0+0.0*pow(gal->MassStar/2.0e10, 2)) + 1.*(Par.SNLoadingFactor-p0) * ( 0.5 *  (1.+gsl_sf_erf((z - 2.5)/1.0)));
				//f = Par.SNLoadingFactor + 3 * ( 0.5 *  (1.+gsl_sf_erf((z - 1.5)/1.0))); // a good model for samples
				//f = Par.SNLoadingFactor + 3 * ( 0.5 *  (1.+gsl_sf_erf((z - 2.5)/1.0)));
				f = Par.SNLoadingFactor + 0 * ( 0.5 *  (1.+gsl_sf_erf((z - 1.2)/1.)));
	/*
	r = dmin(gal->RadiusHalfStar, gal->RadiusHalo);
	sigma0 = gal->RateStarFormation / (2.*M_PI*r*r/1.678/1.678)/1e12;
	if(Do_preheating) f= 1.0*pow(sigma0/300, 6)+Par.SNLoadingFactor * ( 0.5 *  (1.+gsl_sf_erf((z - 2.5)/1.0)));
	*/
	//if (Do_preheating) f = 0.0;
	else f = Par.SNLoadingFactor * pow(200/gal->VelocityVirial, Par.SNLoadingFactorIndex);
	//printf("outflow: %d %g %g %g %g \n", Do_preheating, z, gal->MassStar, r, f);

	return f;
}

double outflow_massloading_factor_surface(double sigma)
{
	double f;

	f = 0 * pow(sigma*1.e-12/2000, 6);

	return f;
}


void star_formation_surface_molecule(struct galaxy *gal, double t, double dt)
{
    int i, j, k, nbin;
    double z0, u0, t0;
    double dden, area, sfe, sfe1, sfr, sd_sfr, ofr, sd_ofr, sd_cold_tmp, mstar, mcold, fml, mp;
    double loadingfactor, tdyn, sfe_max=1e33; //0.4e12;
    double zcold, zyeildhot;

    nbin = gal->nbin;


    u0 = UnitMass_in_g / (UnitLength_in_cm * UnitLength_in_cm);

	j = floor(t/Bin_size_time);
	if (j > N_TIME_BIN) j=N_TIME_BIN;
    tdyn = 1;
    loadingfactor = outflow_massloading_factor(gal);
    for (i=0, sfr=0.0, ofr=0.0, mstar=0.0, mcold=0.0; i<nbin; i++)
    {
    	dden = gal->SDensityColdMolecular[i];
    	//if (gal->MassHalo/mah(0.0) < 1./10) dden = gal->SDensityCold[i];   //**************************************//
    	//if (gal->MassHalo/mah(0.0) < 2./halo_concentration_prada(0, mah(0.0))) dden = gal->SDensityCold[i];
    	sd_cold_tmp = gal->SDensityCold[i];
    	area = M_PI * (gal->RadiusOuter[i]*gal->RadiusOuter[i] - gal->RadiusInner[i] * gal->RadiusInner[i]);

	if ((sd_cold_tmp * area) > 0.0) zcold = (gal->MassMetalCold[i])/(sd_cold_tmp * area);
	else zcold = 0.0;

    	if(dden > 0.0 && sd_cold_tmp >0.0)
    	{
    		//sd_cold_tmp = gal->SDensityCold[i];

    		// model with readjusted star formation efficiency
    		//sfe_max = dden/(gal->SDensityCold[i] * (1+Par.SNLoadingFactor)*dt);
    		// model with readjusted outflow loading factor
    		sfe_max = 1e33;

    		sfe1 = star_formation_efficiency(sd_cold_tmp, tdyn); // depends on total cold gas density in the Krumholz model
    		sfe = dmin(sfe_max, sfe1);
    		//if (gal->MassHalo/mah(0.0) < 1./10) sfe = 1./0.2; //**************************************//
    		//if (gal->MassHalo/mah(0.0) < 2./halo_concentration_prada(0, mah(0.0))) sfe = 1./0.2;
    		sd_sfr = sfe * dden;
    	}
    	else
    	{
    		sfe = sd_sfr = 0.0;
    	}

    		//printf("debug_star_formation: %d %g %g %g %g %g %g\n", i, sd_sfr, sfe1, sfe, dden, sd_cold_tmp);
    		gal->SDensitySFH[i][j] += sd_sfr * dt;
    		for(k=0, t0=0, fml=0, mp=0; k<j; k++)
    		{

    			if(t>gal->TimeArray[0]) fml += gal->SDensitySFH[i][k] *
    					(stellar_mass_loss(t - gal->TimeArray[k]+0.5*dt) -
    							stellar_mass_loss(t - gal->TimeArray[k]-0.5*dt));
				/*
    			if(t>gal->TimeArray[0]) fml += gal->SDensitySFH[i][k] * (stellar_mass_loss(t-gal->TimeArray[k]+0.5*Bin_size_time)-
    							stellar_mass_loss(t-gal->TimeArray[k]-0.5*Bin_size_time));
								*/
    			//t0 = gal->TimeArray[k];
    			mp += gal->SDensitySFH[i][k] * (1.-stellar_mass_loss(t-gal->TimeArray[k]));
    			//if(gal->z < 0.0001) printf("SF debug: %d %d %d %g %g %g %g\n", i, k, j, t, gal->TimeArray[k],stellar_mass_loss(t-gal->TimeArray[k]), mp);
    		}
    		//fml = sd_sfr * dt * 0.43;

    		//printf("%g %g\n", fml, 0.3*sd_sfr*dt);
    		//gal->SDensityStar[i] += (1.-Par.StellarMassLossFraction) * sd_sfr *dt;
    		//gal->SDensityStar[i] -= fml;
    		gal->SDensityStar[i] = mp;

    		sd_cold_tmp -= (1.-Par.StellarMassLossFraction) * sd_sfr *dt;
    		//sd_cold_tmp += fml;

    		//loadingfactor += outflow_massloading_factor_surface(gal->SDensityStar[i]);
    		sd_ofr = dmin(loadingfactor * sd_sfr * dt, sd_cold_tmp);
    		sd_cold_tmp -= sd_ofr;
//printf("debug: %d %g %g %g %g %g\n", i, gal->z, loadingfactor, sd_sfr, dt, fml);
//exit(0);
    		gal->SDensityCold[i] = dmax(0.0, sd_cold_tmp);
    		gal->SDensitySFR[i] = sd_sfr / dt;
			gal->SDensityOFR[i] = sd_ofr / dt;
    		sfr += sd_sfr * area; //if global metallicity calculated in loop, is this necessary?
    		ofr += sd_ofr * area /dt; //if global metallicity calculated in loop, is this necessary?

/*
    		mstar += gal->SDensityStar[i] * area;
    		mcold += gal->SDensityCold[i] * area;
*/
    	//}
    	mstar += gal->SDensityStar[i] * area;
    	mcold += gal->SDensityCold[i] * area;

    // metallicity by radius //sfr and ofr are global so using sd_sfr * area and sd_ofr * area / dt(???) instead
    gal->MassMetalStar[i] += sd_sfr * area * dt * zcold * (1-Par.StellarMassLossFraction);
    gal->MassMetalCold[i] += sd_sfr * area * dt * Par.Yield * (1 - Par.ZFractionYieldToEject - Par.ZFractionYieldToHot) 
      - sd_ofr * area / dt * dt  * zcold - sd_sfr * area * dt * zcold * (1. - Par.StellarMassLossFraction); 
    // gal->MassMetalHot[i] += sdr_sfr * area * dt * Par.Yield * Par.ZFractionYieldToHot; //Assumes hot mass deposited at location of star formation


    //if ((sd_cold_tmp * area) > 0.0) zcoldglobal = (gal->MassMetalCold[i])/(sd_cold_tmp * area);
    //else zcoldglobal = 0.0;

    
    //metallicity global
    gal->MetalStar += sd_sfr * area * dt * zcold * (1.-Par.StellarMassLossFraction);
    gal->MetalCold += sd_sfr * area * dt * Par.Yield * (1. - Par.ZFractionYieldToEject - Par.ZFractionYieldToHot) 
      - sd_ofr * area / dt * dt * zcold - sd_sfr * area * dt * zcold * (1. - Par.StellarMassLossFraction); 
    gal->MetalHot += sd_sfr * area * dt * Par.Yield * Par.ZFractionYieldToHot;
    gal->MetalEject += sd_ofr * area / dt * dt * zcold + sd_sfr *area * dt * Par.Yield * Par.ZFractionYieldToEject;
    //coldadd = sd_sfr * area * dt * Par.Yield * (1 - Par.ZFractionYieldToEject - Par.ZFractionYieldToHot) 
    //  - sd_ofr * area * dt / dt * zcold - sd_sfr * area * dt * zcold * (1. - Par.StellarMassLossFraction);
    //printf("debug_star_formation3: %g %g %g %g\n", t, gal->MetalCold, gal->MassMetalCold[i], coldadd);
    

    }
    //printf("debug_star_formation: %g %g %g %g %g %g\n", t, hubble_time(0) * xH0_recip, dt, sfr, mstar, mcold);
    gal->MassCold = mcold;
    gal->MassStar = mstar;
    gal->MassEject += ofr * dt * (1.-Par.MassFractionEjectToHot);

    gal->RateStarFormation = sfr;
    gal->StarFormationHistory[j] += sfr * dt;
    gal->RateOutflow = ofr;

    


    //printf("SF:%g %g %g %g %g %g\n", gal->MassCold, gal->MassStar, gal->MetalCold, gal->MetalStar, sfr, ofr);
}

double epsilon_disk_guo2011(double vmax)
{
	double e0, vreheat, beta1, edisk;

	//Guo et al. 2013

	e0 = 4;
	vreheat = 80.;
	beta1 = 3.2;

	//Henriques et al. 2013

	e0 = 2.1;
	vreheat = 405;
	beta1 = 0.92;


	edisk = e0 * (0.5 + 1./pow(vmax/vreheat, beta1));
	return edisk;
}

double epsilon_halo_guo2011(double vmax)
{
	double eta0, veject, beta2, ehalo;

	//Guo et al. 2013

	eta0 = 0.18;
	veject = 90.;
	beta2 = 3.2;

	//Henriques et al. 2013

	eta0 = 0.65;
	veject = 336;
	beta2 = 0.46;


	ehalo = eta0 * (0.5 + 1./pow(vmax/veject, beta2));
	return ehalo;
}

double massloading_factor_reheat_guo2011(struct galaxy *gal)
{
	double vmax, edisk, ehalo, elimit, mlf;
	vmax = gal->VelocityMax;
	edisk = epsilon_disk_guo2011(vmax);
	ehalo = epsilon_halo_guo2011(vmax);

	elimit = ehalo * 630.*630 / (vmax * vmax);
	mlf = dmin(edisk, elimit) * ( 1. - 0.43); // 0.43 comes from assumed IMF

	return mlf;
}


double massloading_factor_eject_guo2011(struct galaxy *gal)
{
	double vmax, vvir, edisk, ehalo, elimit, mlf;

	vmax = gal->VelocityMax;
	vvir = gal->VelocityVirial;
	edisk = epsilon_disk_guo2011(vmax);
	ehalo = epsilon_halo_guo2011(vmax);

	elimit = ehalo * 630.*630 / (vvir * vvir);

	mlf = dmax((elimit - edisk), 0.0) * (1. - 0.43); // 0.43 comes from assumed IMF
	return mlf;
}

void star_formation_surface_molecule_with_guo2011_feedback(struct galaxy *gal, double t, double dt)
{
    int i, j, k, nbin;
    double z0, u0, t0;
    double dden, area, sfe, sfe1, sfr, sd_sfr, ofr_reheat, ofr_eject, sd_ofr_reheat, sd_ofr_eject, sd_cold_tmp, mstar, mcold, fml, mp;
    double loadingfactor_reheat, loadingfactor_eject, tdyn, sfe_max=1e33; //0.4e12;
    double zcold;

    nbin = gal->nbin;

    u0 = UnitMass_in_g / (UnitLength_in_cm * UnitLength_in_cm);
    j = floor(t/Bin_size_time);
    if (j > N_TIME_BIN) j=N_TIME_BIN;
    tdyn = 1;
    loadingfactor_reheat = massloading_factor_reheat_guo2011(gal);
    loadingfactor_eject = massloading_factor_eject_guo2011(gal);
    for (i=0, sfr=0.0, ofr_reheat=0.0, ofr_eject=0.0, mstar=0.0, mcold=0.0; i<nbin; i++)
    {

    	dden = gal->SDensityColdMolecular[i];
    	sd_cold_tmp = gal->SDensityCold[i];
    	area = M_PI * (gal->RadiusOuter[i]*gal->RadiusOuter[i] - gal->RadiusInner[i] * gal->RadiusInner[i]);

	if ((sd_cold_tmp * area) > 0.0) zcold = (gal->MassMetalCold[i])/(sd_cold_tmp * area);
	else zcold = 0.0;

    	if(dden > 0.0 && sd_cold_tmp >0.0)
    	{
    		//sd_cold_tmp = gal->SDensityCold[i];

    		// model with readjusted star formation efficiency
    		//sfe_max = dden/(gal->SDensityCold[i] * (1+Par.SNLoadingFactor)*dt);
    		// model with readjusted outflow loading factor
    		sfe_max = 1e33;

    		sfe1 = star_formation_efficiency(sd_cold_tmp, tdyn); // depends on total cold gas density in the Krumholz model
    		sfe = dmin(sfe_max, sfe1);
    		sd_sfr = sfe * dden;
    	}
    	else
    	{
    		sfe = sd_sfr = 0.0;
    	}

    		//printf("debug_star_formation: %d %g %g %g %g %g %g\n", i, sd_sfr, sfe1, sfe, dden, sd_cold_tmp);
    		gal->SDensitySFH[i][j] += sd_sfr * dt;
    		for(k=0, t0=0, fml=0, mp=0; k<j; k++)
    		{

    			if(t>gal->TimeArray[0]) fml += gal->SDensitySFH[i][k] *
    					(stellar_mass_loss(t - gal->TimeArray[k]+0.5*dt) -
    							stellar_mass_loss(t - gal->TimeArray[k]-0.5*dt));
				/*
    			if(t>gal->TimeArray[0]) fml += gal->SDensitySFH[i][k] * (stellar_mass_loss(t-gal->TimeArray[k]+0.5*Bin_size_time)-
    							stellar_mass_loss(t-gal->TimeArray[k]-0.5*Bin_size_time));
								*/
    			//t0 = gal->TimeArray[k];
    			mp += gal->SDensitySFH[i][k] * (1.-stellar_mass_loss(t-gal->TimeArray[k]));
    		}
    		//fml = sd_sfr * dt * 0.43;

    		//printf("%g %g\n", fml, 0.3*sd_sfr*dt);
    		//gal->SDensityStar[i] += (1.-Par.StellarMassLossFraction) * sd_sfr *dt;
    		//gal->SDensityStar[i] -= fml;
    		gal->SDensityStar[i] = mp;

    		sd_cold_tmp -= (1.-Par.StellarMassLossFraction) * sd_sfr *dt;
    		//sd_cold_tmp += fml;
//printf("debug: %d %g %g %g %g %g %g %g\n", i, gal->z, loadingfactor_reheat, loadingfactor_eject, sd_sfr, dt, (loadingfactor_reheat+loadingfactor_eject) * sd_sfr * dt, sd_cold_tmp);

    		//loadingfactor += outflow_massloading_factor_surface(gal->SDensityStar[i]);
    		sd_ofr_reheat = dmin(loadingfactor_reheat * sd_sfr * dt, sd_cold_tmp);
    		sd_cold_tmp -= sd_ofr_reheat;
    		sd_ofr_eject = dmin(loadingfactor_eject * sd_sfr * dt, sd_cold_tmp);
    		sd_cold_tmp -= sd_ofr_eject;

//printf("debug: %d %g %g %g %g %g %g %g\n", i, gal->z, loadingfactor_reheat, loadingfactor_eject, sd_sfr, dt, (loadingfactor_reheat+loadingfactor_eject) * sd_sfr * dt, sd_cold_tmp);
//exit(0);
    		gal->SDensityCold[i] = dmax(0.0, sd_cold_tmp);
    		gal->SDensitySFR[i] = sd_sfr / dt;
			gal->SDensityOFR[i] = (sd_ofr_reheat + sd_ofr_eject) ;
    		sfr += sd_sfr * area; //necessary if global metallicity calculation is in loop?
    		ofr_reheat += sd_ofr_reheat * area /dt;  //necessary if global metallicity calculation is in loop?
    		ofr_eject += sd_ofr_eject * area / dt;  //necessary if global metallicity calculation is in loop?
		//if(gal->SDensityColdMolecular[i] > gal->SDensityCold[i]) printf("%d %g %g\n",i, gal->SDensityCold[i],gal->SDensityColdMolecular[i]);
/*
    		mstar += gal->SDensityStar[i] * area;
    		mcold += gal->SDensityCold[i] * area;
*/
    	//}
    	mstar += gal->SDensityStar[i] * area;
    	mcold += gal->SDensityCold[i] * area;

    //metallicity by radius //sfr, ofr_reheat, and ofr_eject are global so using sd_sfr * area, sd_ofr_reheat * area / dt(???), and  sd_ofr_eject * area / dt(???)
    gal->MassMetalStar[i] += sd_sfr * area * dt * zcold * (1.-Par.StellarMassLossFraction);
    gal->MassMetalCold[i] += sd_sfr * area * dt * Par.Yield * (1. - Par.ZFractionYieldToEject - Par.ZFractionYieldToHot) 
      - (sd_ofr_reheat + sd_ofr_eject) * area / dt * dt * zcold - sd_sfr * area * dt * zcold * (1.-Par.StellarMassLossFraction);
    // gal->MassMetalHot[i] += sd_sfr * area * dt * Par.Yield * Par.ZFractionYieldToHot + sd_ofr_reheat * area * dt * zcold //assumes hot metals deposited locally

    // metallicity global
    //gal->MetalStar += sfr * dt * zcold * (1-Par.StellarMassLossFraction);
    gal->MetalStar += sd_sfr * area * dt * zcold * (1.-Par.StellarMassLossFraction);
    gal->MetalCold += sd_sfr * area * dt * Par.Yield * (1. - Par.ZFractionYieldToEject - Par.ZFractionYieldToHot) 
      - (sd_ofr_reheat + sd_ofr_eject) * area / dt * dt * zcold - sd_sfr * area * dt * zcold * (1.-Par.StellarMassLossFraction);
    gal->MetalHot += sd_sfr * area * dt * Par.Yield * Par.ZFractionYieldToHot + sd_ofr_reheat * area / dt * dt * zcold;
    gal->MetalEject += sd_ofr_eject * area / dt * dt * zcold + sd_sfr * area * dt * Par.Yield * Par.ZFractionYieldToEject;
    gal->MassEject += sd_ofr_eject * area / dt * dt;
    //printf("debug_star_formation: %g\n", zcold);


    }
    //    printf("debug_star_formation: %g %g %g %g\n", t, sfr, mstar, mcold);
    gal->MassCold = mcold;
    gal->MassStar = mstar;
    

    gal->RateStarFormation = sfr;
    gal->RateOutflow = ofr_reheat + ofr_eject;


}

double reincorporation_model_henriques2013(double mvir)
{
	double t0 = 18.0; // in Gyr
	double m0 = 1e10; // in Msun
	double t;

	//t = t0 * m0 / mvir;
	t = Par.ReincorporationTimeScale * m0 / mvir;
	return t;
}

void ejected_gas_reincorporation(struct galaxy *gal, double t, double dt)
{
	double zeject, trein, dm;

	if (gal->MassEject > 1e-33) zeject = gal->MetalEject / gal->MassEject;
	else zeject = 0.0;

	trein = reincorporation_model_henriques2013(gal->MassHalo);

	dm = dmin(gal->MassEject * dt / trein, gal->MassEject);

	gal->MassEject -= dm;
	gal->MetalEject -= dm * zeject;
	gal->MetalHot += dm * zeject;

	//printf("debug_star_formation: %g\n", zeject);
}
