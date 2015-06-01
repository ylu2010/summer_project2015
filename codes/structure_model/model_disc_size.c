/*
 * model_disc_size.c
 *
 *  Created on: Apr 8, 2013
 *      Author: luyu
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

#include "variables.h"
#include "proto.h"
#include "cosmo.h"


double mdisk_exponential(double x)   /* mass of the disk, scaled to halo total mass ; r is scaled to rd */
{
	double mass;

	mass=(1-(1+x)*exp(-x));
	return mass;
}

double fun_fc_disk(double c)
{
	double term1,term2,term3,term;

	term1=(1+c)*(1+c);
	term2=log(1+c);
	term3=(c/(1+c)-term2)*(c/(1+c)-term2);

	term=0.5*c*(1-1/term1-2*term2/(1+c))/term3;
	return term;
}

double fun_fc_disk_approx(double c)
{
	double fc1;

	fc1=2.0/3+pow((c/21.5),0.7);
	return fc1;
}

double rotcurv(double r, double rhalo, double rd, double mhalo, double mdisk, double c)
{
	double rrd, vcDM2, vcDSK2, vc;

	rrd=r/rd;
	vcDM2=G * mhalo * mass_profile_nfw(r/rhalo, c)/r;
	vcDSK2=G * mdisk * mdisk_exponential(rrd)/r;
	vc = sqrt(vcDM2 + vcDSK2);
	return vc;
}

double f_disk(double u, void * params)
{

	double *p;
	double v200, rd, vc, rhalo, mhalo, mdisk, c;
	double term1, term;

	p = (double *)params;
	v200 = p[0]; rd = p[1]; rhalo = p[2]; mhalo = p[3]; mdisk = p[4]; c = p[5];

	term1=exp(-u)*u*u;
	vc=rotcurv(rd*u, rhalo, rd, mhalo, mdisk, c);
	term=term1*vc/v200;
	return term;
}

double fun_fr_disk(double vhalo,double rd, struct galaxy *gal)
{
	double fr;
	double result,error;
	double param[6];

	gsl_integration_workspace *w
		= gsl_integration_workspace_alloc(1000);
	gsl_function F;

	param[0]=vhalo; param[1]=rd; param[2]=gal->RadiusHalo; param[3]=gal->MassHalo; param[4]=gal->MassStar + gal->MassCold; param[5]=gal->ConcenHalo;
	F.function = &f_disk;
	F.params = param;

	//printf("%g %g %g %g\n", gal->RadiusHalfCold , gal->MassCold ,gal->RadiusHalfStar, gal->MassStar);
	//printf("%g %g %g %g %g %g\n", vhalo, rd, param[2], param[3], param[4], param[5]);
	gsl_integration_qagiu(&F, 0, 0, 1e-7, 1000, w, &result, &error);
	fr=2/result;

	gsl_integration_workspace_free(w);
	return fr;
}

double fun_fr_disk_approx(double lam_p, double md, double c)
{
	double term1, term2, term3, f;

	term1=-0.06+2.71*md+0.0047/lam_p;
	term2=1-3*md+5.2*md*md;
	term3=1-0.019*c+0.00025*c*c+0.52/c;
	//printf("t1=%g t2=%g t3=%g", term1, term2, term3);
	f = pow(lam_p/0.1,term1)*term2*term3;
	return f;
}

double disk_radius(struct galaxy *gal)
{
	double md, lam_p, fc, fr, rd;

	if ( Do_preheating)
	{
		fc = fun_fc_disk(gal->ConcenHalo);

		lam_p = gal->SpinHalo * gal->SpinCooling;
		//md = (gal->MassCold + gal->MassStar) / gal->MassHalo;
		//fr = fun_fr_disk_approx(lam_p, md, gal->ConcenHalo); // not self-consistent when computing for newly accreted gas mass
		/*
		if(gal->MassCold+gal->MassStar > 0)
			fr = fun_fr_disk(gal->VelocityVirial, (gal->RadiusHalfCold * gal->MassCold + gal->RadiusHalfStar * gal->MassStar)/(gal->MassCold+gal->MassStar)/1.678, gal);
		else
			fr = 1;
			*/
			fr = 1.0;
		rd = Par.DiskRadiusFactor / sqrt(2*fc) * lam_p * gal->RadiusHalo * fr;// * (1.+Redshift)/2.;
	}
	else
		rd = Par.DiskRadiusFactor * 0.03 * gal->RadiusHalo;

	//printf("disk: %g md=%g lam_p=%g c=%g fc=%g sc=%g\n", Redshift, md, lam_p, gal->ConcenHalo,  fc, gal->SpinCooling);

	return rd;
}


/*

double adiabatic(double ri, void * params)
{
	double r,rd,md,y;
	double *p;
	md= MassDisk;
	p = (double *)params;
	r = p[0]; rd = p[1];

	y=mDSK(r/rd)*r+mHalo_orig(ri)*(1.0-md)*r-mHalo_orig(ri)*ri;

	//printf("in adiabatic r=%g rd=%g ri=%g y=%g\n",r,rd,ri,y);
	//printf("in adiabatic mDSK=%g mhalo=%g\n",mDSK(r/rd),mHalo_orig(ri));

	return y;
}

double mTOT(double r, double rd)
{
	double mass;
	double params[2];
	int status;
	int iter=0, max_iter=1000;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double ri=0;
	double x_low=0.000, x_hig=r*1e6;
	double md=MassDisk;

	gsl_function F;

	params[0] = r;  params[1] = rd;
	F.function = &adiabatic;
	F.params = params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, x_low, x_hig);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		ri = gsl_root_fsolver_root(s);
		x_low = gsl_root_fsolver_x_lower(s);
		x_hig = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_low, x_hig, 0, 0.001);

		//if(status == GSL_SUCCESS) printf("Converged! stu=%d %d\n",GSL_CONTINUE, iter);
	}while(status == GSL_CONTINUE && iter < max_iter);
	mass=mDSK(r/rd)+mHalo_orig(ri)*(1-md);

	gsl_root_fsolver_free(s);
	return mass;
}


double rotcurv(double r, double rd)
{
	double rrd, vcDM2, vcDSK2, vc;

	rrd=r/rd;
	vcDM2=G*(mTOT(r,rd)-mDSK(rrd))/r;
	vcDSK2=G*mDSK(rrd)/r;
	vc = sqrt(vcDM2 + vcDSK2);
	return vc;
}

double f(double u, void * params)
{

	double *p;
	double v200, rd, vc;
	double term1, term;

	p = (double *)params;
	v200 = p[0]; rd = p[1];

	term1=exp(-u)*u*u;
	vc=rotcurv(rd*u,rd);
	term=term1*vc/v200;
	return term;
}

double funcR(double v200,double rd)
{
	double fr;
	double result,error;
	double param[2];

	gsl_integration_workspace *w
		= gsl_integration_workspace_alloc(1000);
	gsl_function F;

	param[0]=v200; param[1]=rd;
	F.function = &f;
	F.params = param;

	gsl_integration_qagiu(&F, 0, 0, 1e-7, 1000, w, &result, &error);
	fr=2/result;

	gsl_integration_workspace_free(w);
	return fr;
}

double funcRfit(double lam, double md, double c)
{
	double term1,term2,term3,term;

	term1=-0.06+2.71*md+0.0047/lam;
	term2=1-3*md+5.2*md*md;
	term3=1-0.019*c+0.00025*c*c+0.52/c;

	term=pow(lam/0.1,term1)*term2*term3;
	return term;
}

double disksize(double jd, double md, double spin, double fc, double c, double vc)
{
	double fr, rd, rd0;
	double lam_p;
	double fr1;

	printf("%f %f %f %f %f %f\n",jd,md,spin,fc,c,vc);
	MassDisk=md;
	Concentration=c;
	fr=1.0;
	rd=1/sqrt(2*fc)*jd/md*spin*fr;
	do
	{
		rd0=rd;

		fr=funcR(vc,rd);

		lam_p=jd/md*spin;
		fr1=funcRfit(lam_p,md,c);
		rd=1/sqrt(2*fc)*jd/md*spin*fr;
		printf("in disksize rd0=%g rd=%g\n",rd0,rd);
	}while(fabs(rd0-rd)/rd >=0.00001);

	printf("in disksize rd=%g fr=%g fr1=%g\n",rd,fr,fr1);
	return rd;
}
*/
