/*
 * model_adiabatic_contraction.c
 *
 *  Created on: Feb 17, 2014
 *      Author: luyu
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "variables.h"
#include "proto.h"
#include "cosmo.h"

struct galaxy *tgal;

double mass_initial(double r)
{
	double m;

	m = mass_profile_nfw(r/tgal->RadiusHalo, tgal->ConcenHalo);
	//m = mass_profile_nfw(r/tgal->RadiusHalo, 5.0);
//	printf("mass_initital: %g %g %g\n", r, tgal->RadiusHalo, m);

	return m;
}

double mass_disk(double r)
{
	double m1, m2, m;

	m1 = interpolate_bipoint(tgal->RadiusOuter, tgal->MassProfCold, N_RADIUS_BIN, r, 0) / tgal->MassHalo;
	m2 = interpolate_bipoint(tgal->RadiusOuter, tgal->MassProfStar, N_RADIUS_BIN, r, 0) / tgal->MassHalo;
	m = m1 + m2;

	return m;
}

double mass_final(double r)
{
	double m;

	m = mass_initial(r) * (1.-BaryonFrac) + mass_disk(r);
	//printf("mass_final: %g %g %g\n", r, mass_initial(r), mass_disk(r));

	return m;
}

double adiabatic_contraction_eqn(double logr1, void *params)
{
    double logr2 = *(double *)params;
    double m1, m2, r1, r2;

    r2 = exp(logr2);
    m2 = mass_final(r2);
    r1 = exp(logr1);
    m1 = mass_initial(r1);
//printf("adiabatic_contraction_eqn: %g %g %g %g %g %g\n", logr1+log(m1), logr2+log(m2), logr1, logr2, log(m1), log(m2));
    return logr1+log(m1) - logr2 - log(m2);
}

void adiabatic_contraction(struct galaxy *gal)
{
	int status,i;
	int iter = 0, max_iter = 10000;
	const gsl_root_fsolver_type *T;
        gsl_root_fsolver *s;
    double logr;
    double x_low,x_hig;
    double root, r_initial;
    tgal = gal;

    gsl_function F;
        F.function = &adiabatic_contraction_eqn;

    for(i=0;i<gal->nbin;i++)
    {
    	logr = log(gal->RadiusOuter[i]);
    	F.params = &logr;

    	iter = 0;
        x_low = log(gal->RadiusMostInner)-1; x_hig = log(gal->RadiusMostOuter);
        T = gsl_root_fsolver_brent;
        s = gsl_root_fsolver_alloc(T);
        gsl_root_fsolver_set (s, &F, x_low, x_hig);

        do{
        	iter++;
            status = gsl_root_fsolver_iterate(s);
            root = gsl_root_fsolver_root(s);
            x_low = gsl_root_fsolver_x_lower(s);
            x_hig = gsl_root_fsolver_x_upper(s);
// printf("low=%g hig=%g\n",x_low,x_hig);
            status = gsl_root_test_interval(x_low, x_hig, 0, 1.e-4);
        }while (status == GSL_CONTINUE && iter < max_iter);
        r_initial = exp(root);
        //printf("adiabatic_contraction: %g %g %g %g\n", logr, root, gal->RadiusOuter[i], r_initial);
        gal->MassProfDMContracted[i]= gal->MassHalo * (1. - BaryonFrac) * mass_initial(r_initial);
        gal->MassProfDM[i] = gal->MassHalo * (1. - BaryonFrac) * mass_initial(gal->RadiusOuter[i]);
//        printf("adiabatic_contraction: %g %g %g %g %g %g\n", logr, root, gal->RadiusOuter[i], r_initial, gal->MassProfDMContracted[i], gal->MassProfDM[i]);
        gsl_root_fsolver_free(s);
    }
}
