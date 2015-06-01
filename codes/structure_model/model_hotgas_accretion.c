/*
 * model_hotgas_accretion.c
 *
 *  Created on: Apr 22, 2013
 *      Author: luyu
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "variables.h"
#include "proto.h"
#include "cosmo.h"


double preheating(double z, double zph, double dzph)
{
    double y;
    y = 1.- 0.5 *  (1.+gsl_sf_erf((z - zph)/dzph));
    return y;
}



double preheating_halo_mass(double z)
{
	double mc, lgmc, m0=1e12; //m0 = 3e11;
	double p[] = {12.6403,		     -1.81156,		     0.122902,		  -0.00380529};

	double beta = 0.6; //0.4;
	double gamma = 0.9;
	mc = 1e12*pow(1+z, beta) / exp(gamma*z);

	lgmc = log10(mc);

/*
	lgmc = p[0] + p[1] * z + p[2] * z *z + p[3] * z * z * z;
	mc = pow(10., lgmc)/xhubble;
*/
	//printf("preh: z=%g lgmc=%g mc=%g\n",  z, lgmc, mc);
	return mc;
}

double do_preheating(struct galaxy *gal, double zcurr)
{
	double mh_c, rho, v, temp, s_vir, s_ph, ss, fraction, s0;
	double factor;

	//factor = 1.e3 * ELECTRONVOLT * pow(BaryonFrac/0.6/PROTONMASS*UnitMass_in_g, 2./3)/BOLTZMANN / (UnitLength_in_cm * UnitLength_in_cm);

	//mh_c = preheating_halo_mass(zcurr);
	//rho = Delta_vir(zcurr) * rho_crit(zcurr) * xhubble * xhubble;
	//v = sqrt(G) * pow(4./3*M_PI*rho, 1./6) * pow(mh_c, 1./3);
	//v = pow(G*mh_c*xH(z) * sqrt(0.5*Delta_vir(z)), 1./3);  // just another way
	//temp = 35.9 * v * v;
	//s_ph = 0.8*temp/pow(rho, 2./3) / factor;

	//s_ph = 14./pow(1+zcurr, 1.);
	//s_ph = 10.0;
	//if(zcurr > 4) s_ph=0.001;
	//s_ph = 15 * preheating(zcurr, 2, 0.4);
	//s_ph = s_ph / pow(1.+zcurr, 2);
	//s_ph = 15./pow(1+zcurr, 0.);
	//printf("%g %g\n", zcurr, s_ph);
	//s_ph = 10*exp(-zcurr/1.3)*pow(1.+zcurr, 0.3);
	//s_ph = 12*exp(-zcurr/1.5)*pow(1.+zcurr, 0.2);

	//s0 = 12;
	//s0 = 11* pow(mah(0.0)/1e11, 0.2);
	//s0 = 17* pow(mah(0.0)/1e12, 0.2); // this is the model in submitted version !!!!!
	s0 = Par.PreheatEntropy*pow(mah(0.0)/1.e12, Par.PreheatEntropySlope);
	//s_ph = 10*exp(-zcurr/1.3)*pow(1.+zcurr, 0.3);
	//s0 = 50* pow(mah(0.0)/1.e12, 0.3);
	//s0 = 15;
	//s_ph = s0/(1.+pow(zcurr/1.1, 2));
	//s_ph = s0/(1.+pow(zcurr/1.2, 2)); // The disk paper used this
	s_ph = s0/(1.+pow(zcurr/1.2, 2.));
	//s_ph = s0/(1.+pow(zcurr/(1.2*(1.-0.3*log10(mah(0.0)/1.e11))), 2));


	s_vir = gal->EntropyVirial;
    ss = s_ph/s_vir/0.8;
    //Par.EntropyRatio = (s_ph/s_vir);
    //Par.EntropyRatio = dmax(s_ph/s_vir, 1.0);

    if(s_ph < s_vir)
    {
    	//Par.EntropyProfileIndex = 1.1;  // Assuming beta===0 in preheating scenario no matter if S_vir < S_ph or not.
    	Par.EntropyRatio = 1.0;
    }
    else
    {
    	//Par.EntropyProfileIndex = 0.0;  // Assuming beta===0 in preheating scenario no matter if S_vir < S_ph or not.
    	Par.EntropyRatio = s_ph/s_vir;
    }

    //Par.EntropyRatio = s_ph/s_vir;
    //Par.EntropyRatio  = 1.0;


    fraction = 1./sqrt(1. + ss*ss*ss);
    //fraction = 1.;
//printf("entropy: %g %g f=%g %g %g %g\n", Par.EntropyProfileIndex, Par.EntropyRatio, fraction, s_ph, s_vir, zcurr);
    return fraction;
}

double hot_accretion_fraction(struct galaxy *gal, double z)
{
	double f;
	if(Do_preheating) f = do_preheating(gal, z);
	else f = 1;

	return f;
}


double hot_gas_accretion_fraction(struct galaxy *gal, double z)
{
	double rate;
	double fhot=1;

	fhot = hot_accretion_fraction(gal, z);
	rate = BaryonFrac * Par.BaryonAccretionFraction * do_reionization(gal->MassHalo, z) * fhot;

	//printf("hot accretion: %g %g %g %g\n", do_reionization(gal->MassHalo, z), do_preheating(gal, z), rate, fhot);
	//rate = BaryonFrac * 0.5;
	return rate;
}

double halo_cold_gas_accretion_fraction(struct galaxy *gal, double z)
{
	double frac=0.0;
	double fhot;
/*
	fhot = hot_accretion_fraction(z, gal->MassHalo);
	frac = BaryonFrac * Par.BaryonAccretionFraction * do_reionization(gal->MassHalo, z) * do_preheating(gal, z) * (1-fhot);
*/
	return frac;
}


double gas_density_profile_r2(double x, double c)
{
	double rho = 1./x/x;
	return rho;
}

double gas_density_profile_maller(double x, double c)
{
	double rho;
	rho = pow(1.+3.7/(x*c)*log(1.+x*c)-3.7/c*log(1+c), 1.5);
	return rho;
}


double gas_temperature_profile_maller(double x, double c)
{
	double t = 1+ 3.7/(x*c)*log(1.+x*c) - 3.7/c*log(1+c);
	return t;
}



double gas_density_profile(double x, double c)
{
	double f;
	//f = gas_density_profile_maller(x, c);
	f = gas_density_profile_r2(x, c);

	return f;
}


double gas_temperature_profile(double x, double c)
{
	double t=1;

	//double t = tvir * gas_temperature_profile_maller(x, c);

	return t;
}

double int_den(double x, void *params)
{
	double *p = (double *) params;
	double c, f;

	c = p[0];

	f = gas_density_profile(x, c) * x * x;
	return f;
}

void hot_gas_profile(struct galaxy *gal)
{
	int i;
	double mvir, rvir, hotgas, vvir, t_gas, den0, den, c, r, x;
	double result, error;
	double p[1];
	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
	gsl_function F;
	F.function = &int_den;
	F.params = p;

	p[0] = c = gal->ConcenHalo;

	mvir = gal->MassHalo;
	rvir = gal->RadiusHalo;
	hotgas = gal->MassHot;
	vvir = gal->VelocityVirial;
	t_gas = 35.9 * vvir*vvir;

	gsl_integration_qags(&F, 0, 1, 0, 1.e-7, 1000, w, &result, &error);
	den0 = hotgas/(4*M_PI*rvir*rvir*rvir*result);

	for(i=0; i<gal->nbin; i++)
	{
		r = 0.5*(gal->RadiusInner[i]+gal->RadiusOuter[i]);
		if(r < gal->RadiusHalo)
		{
			x = r/gal->RadiusHalo;
			gal->DensityProfHot[i] = den0 * gas_density_profile( x, p[0] );

			gsl_integration_qags(&F, 0, x, 0, 1.e-7, 1000, w, &result, &error);
			gal->MassProfHot[i] = 4.*M_PI*(rvir*rvir*rvir)*den0*result;
			//printf("int: %g %g %g %g %g\n", r, x, gal->RadiusHalo, gal->DensityProfHot[i], gal->MassProfHot[i]);
			gal->TemperatureProfHot[i] = t_gas * gas_temperature_profile(x, c);
		}
		else
		{

		}
	}

	gsl_integration_workspace_free(w);
}

double gas_density_profile_isenthermal_model1(double x, double z, double s_ratio)
{
	double delta_c, delta_m, rho_c, rho_c0, zp1, rho0, term1, rho;
/*
	delta_c = Delta_vir(z);
	rho_c = rho_crit(z);
	rho_c0 = rho_crit(0);
	zp1 = z+1;

	delta_m = delta_c * rho_c/rho_c0/Omega_M0/zp1/zp1/zp1;
*/
	term1 = 0.8 * s_ratio * log(x);
	rho = pow(1 - term1, 1.5);
	return rho;
}



double gas_density_profile_isenthermal_model2(double x, double z, double s_ratio, double c)
{
	double delta_c, delta_m, rho_c, rho_c0, zp1, rho0, term1, rho;
	double cx, term_c, term_nfw;
/*
	delta_c = Delta_vir(z);
	rho_c = rho_crit(z);
	rho_c0 = rho_crit(0);
	zp1 = z+1;
	rho0 = rho_c0 * Omega_M0 * zp1*zp1*zp1;// * BaryonFrac;
	delta_m = delta_c * rho_c/rho0;
*/
	cx = c*x;
	term_nfw = log(1+cx)/cx - log(1+c)/c;//- 2*(log(cx/(1+cx)) - log(c/(1+c)));
	term_c = c/(log(1+c)-c/(1+c));

	term1 = 0.8 * s_ratio;

	rho = pow( (1. + term1 * term_c * term_nfw), 1.5);
	return rho;
}

double gas_density_profile_isenthermal(double x, double z, double s_ratio, double c)
{
	double rho;
	rho = gas_density_profile_isenthermal_model1(x, z, s_ratio);
	//rho = gas_density_profile_isenthermal_model2(x, z, s_ratio, c);

	return rho;
}

double gas_temperature_profile_isenthermal(double rho_ratio, double z, double s_ratio, double c)
{
	double t;

	t = pow( rho_ratio, 2./3) / s_ratio;
	return t;
}

double int_den_isenthermal(double x, void * params)
{
	double *p = (double *) params;
	double z, s_ratio, c, f;

	z = p[0];
	s_ratio = p[1];
	c = p[2];

	f = gas_density_profile_isenthermal(x, z, s_ratio, c) * x * x;
	return f;
}

double temperature_igm(double z)
{
	double t;// = pow(10., 0.26+0.92*z-0.23*z*z) - 1;
	t = 72e5/((1.+z)*(1.+z));

	return t;
}

void hot_gas_profile_isenthermal(struct galaxy *gal)
{
	int i;
	double z, c, s_ratio;
	double r, x, mvir, rvir, vvir, t_gas, hotgas;
	double delta_c, rho_c, rho_c0, zp1, rho0, delta_m;
	double t_igm, rhomean, rhovir;
	double result, error;
	double p[3];

	gsl_integration_workspace * w = gsl_integration_workspace_alloc(1000);
	gsl_function F;
	F.function = &int_den_isenthermal;
	F.params = p;
	p[0] = z = gal->z;
	p[2] = c = gal->ConcenHalo;

	mvir = gal->MassHalo;
	rvir = gal->RadiusHalo;
	vvir = gal->VelocityVirial;
	t_gas = 35.9 * vvir*vvir;
	s_ratio = 1./Par.EntropyRatio;
	/*
	if(z<2)
	{

	}
	*/
	p[1] = s_ratio;

	delta_c = Delta_vir(z);
	rho_c = rho_crit(z);
	rho_c0 = rho_crit(0);
	zp1 = z+1;
	rho0 = rho_c0 * Omega_M0 * zp1*zp1*zp1;// * BaryonFrac;
	delta_m = delta_c * rho_c/rho0;

	gsl_integration_qags(&F, 0, 1, 0, 1.e-8, 1000, w, &result, &error);
	//printf("isenthermal: z=%g s_ratio=%g c=%g result=%g\n", p[0], p[1], p[2], 4*M_PI*(rvir*rvir*rvir)*result);
	rhovir = gal->MassHot / (4*M_PI*(rvir*rvir*rvir)*result);

	//printf("isenthermal check: %g %g\n", gal->MassHot, result);

	for(i=0; i<gal->nbin; i++)
	{
		r = 0.5*(gal->RadiusInner[i]+gal->RadiusOuter[i]);
		if(r < gal->RadiusHalo)
		{
			x = r/gal->RadiusHalo;
			gal->DensityProfHot[i] = rhovir * gas_density_profile_isenthermal( x, p[0], p[1], p[2] );

			gsl_integration_qags(&F, 0, x, 0, 1.e-7, 1000, w, &result, &error);
			gal->MassProfHot[i] =  4.*M_PI*(rvir*rvir*rvir)*rhovir*result;

			gal->TemperatureProfHot[i] = t_gas * gas_temperature_profile_isenthermal(gal->DensityProfHot[i]/rhovir, z, s_ratio, c);
			//printf("%g %g %g %g %g %g\n", r, x, gal->DensityProfHot[i]/(rhovir * gas_density_profile_isenthermal( 1, p[0], p[1], p[2] )), gal->MassProfHot[i], gal->TemperatureProfHot[i]);
	    	printf("%g %g %g %g %g %g\n", r, x, gal->RadiusHalo, gal->DensityProfHot[i], gal->MassProfHot[i], gal->TemperatureProfHot[i]);
		}
		else
		{

		}
	}

	gsl_integration_workspace_free(w);
	exit(0);
}

double gas_temperature_profile_power_law_entropy(double x, double rho_ratio, double s_ratio, double beta)
{
	double t;

	t = pow(x, beta) * pow( rho_ratio, 2./3) / s_ratio;
	return t;
}


double m_halo(double x, double c)
{
	double m0, mr;
	// NFW

	m0 = 1.0/(log(1.+c) - c/(1.+c));
	mr = m0 *(log(1.+x*c) - x*c/(1.+x*c));

	//isothermal
	//mr = x;
	return mr;
}

double dmdx(double x, double c)
{
	double m0, dmdx;
	// NFW

	m0 = 1.0/(log(1.+c) - c/(1.+c));
	dmdx = (c/(1+x*c) - ( c/ (1+x*c) - x * c *c /(1+x*c)/(1+x*c))) * m0;

	// isothermal
	//dmdx = 1;
	return dmdx;
}

int func_power_law_entropy(double xp, const double y[], double f[], void *params)
{
	double *p = (double *)params;
	double x, adia, beta, c;
	x = xp;
	adia = p[0]; beta = p[1]; c = p[2];
	//printf("a=%g b=%g\n", adia, beta);
	f[0] = - 0.4 * ( beta / x * y[0] + 2. * adia * m_halo(x, c)/ pow(x, beta+2));
	return GSL_SUCCESS;
}

int jac_power_law_entropy (double xp, const double y[], double *dfdy,  double dfdt[], void *params)
{
	double *p = (double *)params;
	double x, adia, beta, c;

	x = xp;
	adia = p[0]; beta = p[1]; c = p[2];
	gsl_matrix_view dfdy_mat
    	= gsl_matrix_view_array (dfdy, 1, 1);
	gsl_matrix * m = &dfdy_mat.matrix;
	gsl_matrix_set (m, 0, 0, -0.4 * beta / x);
	dfdt[0] = 0.4 * (beta * y[0]/x/x - 2 * adia * (dmdx(x, c) / pow(x, beta+2) - (beta+2) * m_halo(x, c) / pow(x, beta+3)));
	return GSL_SUCCESS;
}

void hot_gas_profile_power_law_entropy ( struct galaxy *gal )
{
	int i;
	double r, x, rho, rhovir, mass, vvir, t_gas, tnorm;
	double s_ratio, beta;

	s_ratio = 1./Par.EntropyRatio;
	beta = Par.EntropyProfileIndex;

	vvir = gal->VelocityVirial;
	t_gas = 35.9 * vvir*vvir;
	//double factor = 1.e3 * ELECTRONVOLT * pow(BaryonFrac/0.6/PROTONMASS*UnitMass_in_g, 2./3)/BOLTZMANN / (UnitLength_in_cm * UnitLength_in_cm);
	double factor = 1.e3 * ELECTRONVOLT * pow(1./0.6/PROTONMASS*UnitMass_in_g, 2./3)/BOLTZMANN / (UnitLength_in_cm * UnitLength_in_cm);

	const gsl_odeiv_step_type * T
			= gsl_odeiv_step_rk8pd;
			//= gsl_odeiv_step_rk4;
			//= gsl_odeiv_step_rk4imp;
			//= gsl_odeiv_step_bsimp;

	gsl_odeiv_step * s
         = gsl_odeiv_step_alloc (T, 1);
	gsl_odeiv_control * c
         = gsl_odeiv_control_y_new (1e-6, 0.0);
	gsl_odeiv_evolve * e
         = gsl_odeiv_evolve_alloc (1);
	int status;

	double params[3] = {s_ratio, beta, gal->ConcenHalo};
    gsl_odeiv_system sys = {func_power_law_entropy, jac_power_law_entropy, 1, &params};

    double x1 = 1.0, x0 = 0, xp;
    double h = -1e-6;
    double y[1] = { 1.0 };

    xp = 1;
    for (i = gal->nbin-1; i >= 0; i--)
    {
    	r = 0.5*(gal->RadiusInner[i]+gal->RadiusOuter[i]);
    	if(r < gal->RadiusHalo)
    	{
    		x = r/gal->RadiusHalo;
    		while ( xp > x )
    		{
    			status = gsl_odeiv_evolve_apply(e, c, s, &sys, &xp, x, &h, y);
    			if(status != GSL_SUCCESS) break;
    		}
    		rho = pow(y[0], 1.5);

    		gal->DensityProfHot[i] = rho;
    		gal->TemperatureProfHot[i] = gas_temperature_profile_power_law_entropy(x, rho, s_ratio, beta);
    		//printf("%d %g %g %g\n", i, x, rho, gal->TemperatureProfHot[i]);
    	}
    	else
    	{
    		gal->DensityProfHot[i] = 0.0;
    		gal->TemperatureProfHot[i] = 0.0;
    	}
    }
    x = gal->RadiusInner[0] / gal->RadiusHalo;
    status = gsl_odeiv_evolve_apply(e, c, s, &sys, &xp, x, &h, y);
    rho = pow(y[0], 1.5);
    //printf("%d %g %g %g\n", i, x, rho, h);
    mass = 4 * M_PI / 3 * rho * gal->RadiusInner[0] * gal->RadiusInner[0] * gal->RadiusInner[0];
    for(i=0; i<gal->nbin; i++)
    {
    	r = 0.5*(gal->RadiusInner[i]+gal->RadiusOuter[i]);
    	mass += 4.*M_PI * gal->DensityProfHot[i] * r * r * (gal->RadiusOuter[i] - gal->RadiusInner[i]);
    	gal->MassProfHot[i] = mass;
    }
    //gal->MassHot = dmax(gal->MassHalo * BaryonFrac * dmin(gal->MassHalo/mah_mcbride(gal->z), 1.) - (gal->MassStar+gal->MassCold), 0);
    //gal->MassHot = dmax(gal->MassHalo * BaryonFrac * dmin(gal->MassHalo/mah_wechsler(1./(gal->z+1)), 0.8) - (gal->MassStar+gal->MassCold), 0);
    //gal->MassHot = dmax(gal->MassHalo * BaryonFrac - (gal->MassStar+gal->MassCold), 0);
    // renormalize the profile according to total hot gas mass
    //if(Do_preheating ) tnorm = t_gas * pow(1*gal->MassHot/mass, 2./3) * gal->EntropyVirial / gal->TemperatureVirial * factor; // when x=1 rho=1
    if(Do_preheating ) tnorm = pow(1*gal->MassHot/mass, 2./3) * gal->EntropyVirial * factor; // when x=1 rho=1
    else tnorm=t_gas;
    //tnorm=t_gas;
    //tnorm = t_gas * pow(1*gal->MassHot/mass, 2./3) * gal->EntropyVirial / gal->TemperatureVirial * factor; // when x=1 rho=1
    for(i=0; i<gal->nbin; i++)
    {
    	r = 0.5*(gal->RadiusInner[i]+gal->RadiusOuter[i]);
    	gal->DensityProfHot[i] *= gal->MassHot/mass;
    	gal->MassProfHot[i] *= gal->MassHot/mass;
    	gal->TemperatureProfHot[i] *= tnorm;


    	//printf("hot_gas: %g %g %g %g %g %g\n", r, x, gal->RadiusHalo, gal->DensityProfHot[i], gal->MassProfHot[i], gal->TemperatureProfHot[i]);
    }

//printf("test: %g %g %g\n", pow(rho1*gal->MassHot/mass, 2./3) ,1/( gal->EntropyVirial / gal->TemperatureVirial), factor);
    //printf("rho1=%g\n", rho1);
    gsl_odeiv_evolve_free (e);
    gsl_odeiv_control_free (c);
    gsl_odeiv_step_free (s);
    //exit(0);
}
