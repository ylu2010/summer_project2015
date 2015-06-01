#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_integration.h>

#include "cosmo.h"


double xH0;
double xH0_recip;
double htime_z0;
double H0;
double xhubble;
double Omega_M0;
double Omega_Lambda0;
double Omega_Baryon0_h2;
double Omega_Baryon0;
double MF_WDM;
double G;
double Sigma8;
double nspec;
int BBSK;
int EBW;
int EISHU;
int ifilter;
int ivar;
double Gamma_cosmo;
double rho_crit0;
double rho_aver0;
double sigma8_norm;

double rfilter, f_bar;
double sEH, bnode, ksilk, keq, alpha_c, alpha_b, beta_c, beta_b;


void init_cosmo(void)
{
	int iftemp, ivtemp;
	double m_WDM, c8, sTH, sGS, sSK;
    double f, bb1, bb2, zeq, zd, y, Fy, Gy, Req, Rd, a1, a2, b1, b2;

	H0 = 100.0;
	xhubble = 0.7;
	Omega_M0 = 0.27;
	Omega_Lambda0 = 0.73;
	Omega_Baryon0_h2 = 0.02233;
	Omega_Baryon0 = 0.041 ;
	MF_WDM = 0;
	G = 4.2994e-9; // [Mpc km^2/Msun/s^2]
	Sigma8=0.79;
	nspec=0.960;
	BBSK=1;
	EBW=0;
	EISHU=0;
	ifilter=1;
	ivar=1;


	xH0 = xhubble * H0;
	xH0_recip = 1.0 / (xH0 * YEAR_IN_SECOND/PC_IN_CM*1.e8);
	htime_z0 = hubble_time(0);
    rho_crit0 = 3.0*H0*H0 / (8.0 * PI * G);
    rho_aver0 = rho_crit0 * Omega_M0;
	Omega_Baryon0 = Omega_Baryon0_h2/xhubble/xhubble;
	f_bar = Omega_Baryon0/Omega_M0;
/*
c---define Gamma: two options here; i) go with the standard definition
c   (Gamma = Omega_0 h), which is for negligible baryon mass, or use
c   the baryonic correction from Sugiyama (1995).

c---with baryonic correction

      Gamma_cosmo = omega_0 * xhubble *
    &               XEXP(-omega_b - SQRT(2.0*xhubble) * f_bar)

c---without baryonic correction

c     Gamma_cosmo = omega_0 * xhubble
c      Gamma_cosmo = 0.21
*/
	//Gamma_cosmo = Omega_M0 * xhubble;
	Gamma_cosmo=0.20;

// if we are considering WDM, compute the mass of WDM particles
    if(MF_WDM!=0.0)
    {
		m_WDM = 2.4 * pow(xhubble, 1.25) * sqrt(Omega_M0-Omega_Baryon0) *
                pow(MF_WDM/1.0e+11, -0.25);
        rfilter = 0.065 * pow(Omega_M0-Omega_Baryon0, -1.0/3.0) *
                  pow(MF_WDM/1.0e+11, 1.0/3.0);
        printf("m_WDM = %g keV\n", m_WDM);
        printf("R_f = %g Mpc/h\n", rfilter);
    }

/*
c---define a number of parameters needed to compute the Eisenstein & Hu
c   power specrum
*/
    f = Omega_M0 * xhubble*xhubble;

    bb1 = 0.313 * pow(f, -0.419) * (1.0 + 0.607*pow(f, 0.674));
    bb2 = 0.238 * pow(f, 0.223);

    bnode = 8.41* pow(f, 0.435);

    keq = 7.46e-2 * f;
    ksilk = 1.6 * pow(Omega_Baryon0_h2, 0.52) * pow(f, 0.73) *
            (1.0 + pow(10.4*f,-0.95));

    zeq = 2.5e+4 * f;
    zd = 1291.0 * pow(f, 0.251)/(1.0 + 0.659*pow(f, 0.828)) *
         (1.0 + bb1 * pow(Omega_Baryon0_h2, bb2));
    y = ((1.0+zeq)/(1.0+zd));
    Fy = log((sqrt(1.0+y) + 1.0) / (sqrt(1.0+y) - 1.0));
    Gy = y * (-6.0 * sqrt(1.0+y) + (2.0+3.0*y) * Fy);

    Req = 31.5 * Omega_Baryon0_h2 * (1000.0/zeq);
    Rd  = 31.5 * Omega_Baryon0_h2 * (1000.0/zd);

    sEH = (2.0/(3.0*keq)) * sqrt(6.0/Req) *
          log((sqrt(1.0+Rd) + sqrt(Rd+Req))/(1.0+sqrt(Req)));

    a1 = pow(46.9*f, 0.670) * (1.0 + pow(32.1*f, -0.532));
    a2 = pow(12.0*f, 0.424) * (1.0 + pow(45.0*f, -0.582));
    b1 = 0.944 / (1.0+pow(458.0*f, -0.708));
    b2 = pow(0.395*f, -0.0266);

    alpha_c = pow(a1, -f_bar) * pow(a2, -(f_bar*f_bar*f_bar));
    beta_c = 1.0 + b1*(pow((Omega_M0-Omega_Baryon0)/Omega_M0, b2) - 1.0);
    beta_c = 1.0/beta_c;

    alpha_b = 2.07 * keq * sEH * pow(1.0+Rd, -0.75) * Gy;
    beta_b = 0.5 + f_bar +
             (3.0-2.0*f_bar) * sqrt(pow(17.2*f, 2) + 1.0);

/*
c---compute un-normalized rms mass variance inside spherical
c   shell with radius Rf = 8 h^{-1} Msun (i.e., filter of choice)
c   This is used to normalize the power-spectrum to sigma8.
c   The mass inside this sphere depends on the filter.
*/
	iftemp = ifilter;
    ivtemp = ivar;
    ivar = 1;
    sigma8_norm = Sigma8;

        switch (ifilter)
        {
                case 1:
                        sigma8_norm = var_numerical(5.9543E+14 * Omega_M0);
                        c8 = 1.0;
                        break;
                case 2:
                        ifilter=1;
                        sigma8_norm = Sigma8;
                        sTH = var_numerical(5.9543E+14 * Omega_M0);
                        ifilter=2;
                        sigma8_norm = Sigma8;
                        sGS = var_numerical(2.2388E+15 * Omega_M0);
                        c8 = sGS/sTH;
                        ifilter=2;
                        sigma8_norm = Sigma8;
                        sigma8_norm = var_numerical(2.2388E+15 * Omega_M0);
                        break;
                case 3:
                        ifilter=1;
                        sigma8_norm = Sigma8;
                        sTH = var_numerical(5.9543E+14 * Omega_M0);
                        ifilter=3;
                        sigma8_norm = Sigma8;
                        sSK = var_numerical(8.4177E+15 * Omega_M0);
                        c8 = sSK/sTH;
                        ifilter=3;
                        sigma8_norm = Sigma8;
                        sigma8_norm = var_numerical(8.4177E+15 * Omega_M0);
                        break;
                default:
                        printf("give an appropriate filter number!\n");
                        exit(2);
        }

        if(ifilter!=iftemp)
        {
                printf("filter stalled up! exit at init_cosmo!\n");
                printf("%d %d\n",ifilter,iftemp);
                exit(2);
        }

        ivar = ivtemp;
        //printf("c8 = %g\n",c8);
}

double xH(double z) 
// calculate hubble constant ( in physical units: km/s/Mpc ) at redshift z.
{
        double z1, fac, xH;

        z1 = 1.0 + z;
        fac = Omega_Lambda0 + (1.0 - Omega_Lambda0 - Omega_M0)*z1*z1 +
                Omega_M0*z1*z1*z1;
        xH = xH0 * sqrt(fac);

        return xH;
}

double Omega_m(double z) 
// calculates the density parameter omega_m at redshift z. 
{
        double omega;

        omega = Omega_M0 * (1.0+z)*(1.0+z)*(1.0+z) / pow(xH(z)/xH0, 2);
        return omega;
}

double rho_crit(double z) // critical density of the universe
{
        double h,rho_c;

        h = xH(z)/xH0 * H0;
        rho_c = 3*h*h / (8*M_PI*G);   //in unit of Msun h^-1/(Mpc h^-1)^3
        return rho_c;
}


double lin_growth(double z) 
// calculates the linear growth factor at redshift z relative to z=0.
// based on Mo's lecture notes and Carroll etal 1992
{

        double o_m, o_l, g, g0, D;

        o_m = Omega_m(z);
        o_l = Omega_Lambda0 / (xH(z)/xH0 * xH(z)/xH0);

        g = 2.5 * o_m / (pow(o_m, 4.0/7.0) - o_l + (1.0 + 0.5*o_m)*(1.0 + o_l/70.0));
	o_m = Omega_M0;
	o_l = Omega_Lambda0;
	g0 = 2.5 * o_m / (pow(o_m, 4.0/7.0) - o_l + (1.0 + 0.5*o_m)*(1.0 + o_l/70.0));
        D = g/g0/(1 + z);

        return D;
}


double delta_c(double z)
// Calculate the critical overdensity.
// We use the approximation in NFW97
// copied from Frank's program
{
	double ddd,dc0,omz;
      	double growth_rate, delta_c;

      	omz = Omega_m(z);
      	ddd = 1.0 / lin_growth(z);
      	dc0 = 0.0;

      	if (fabs(1.0-Omega_M0-Omega_Lambda0) < 1.0E-4)
        	dc0 = 0.15 * pow(12.0 * PI, 2.0/3.0) * pow(omz, 0.0055);

      	if (Omega_M0 < 1.0 && Omega_Lambda0 == 0.0) 
        	dc0 = 0.15 * pow(12.0 * PI, 2.0/3.0) * pow(omz, 0.0185);

      	if (dc0 == 0.0)
	{
        	printf("delta_c is not defined for this cosmology! \n");
       		printf(" omega_0, omega_lambda: %g %g\n",Omega_M0,Omega_Lambda0);
        	printf(" omega_z, growth_rate: %g %g\n", omz,ddd);
        	printf(" at redshift: %g\n", z);
		exit(10);
	}
      	else
	{
        	delta_c = dc0 * ddd;
	}
	return delta_c;
}

double Delta_vir(double z)
// calculate the virial density in terms of critical density of the Universe.
// We use the fitting formulae of Bryan & Norman, 1998, ApJ, 495, 80.
// These are accurate to better than 1% for 0.1<Omega_M0<1.0
{
	int flag = 0;
	double x, dc;

	x = Omega_m(z) -1.0;
        if(Omega_M0 < 1.0 && Omega_Lambda0==0.0)
	{
                dc = 18.0*PI*PI + 60.0*x - 32.0*x*x;
		flag = 1;
	}
       	if(fabs(1.0-Omega_M0-Omega_Lambda0) < 1.e-4)
	{
               	dc = 18.0*PI*PI + 82.0*x - 39.0*x*x;
		flag = 2;
	}
	if(!flag)
	{
		printf("Delta_cric(z) is not defined for this cosmology!\n Program exits!\n");
		exit(10);
	}
	return dc;
}


double hubble_time(double z)
// calculate hubble time at redshift z in units of (1/H0).
// It has been check for SCDM, OCDM, and LCDM, where it works fine, by Frank
{
	int flag2;
	double tz, z1, a, atmp, omega1, const0, param0;
	double xx, xxz, xxmax, step;

	if(Omega_M0 == 1.0)
	{
		tz = (2.0/3.0)*pow(1.+z, -1.5);
	//printf("f=1\n");
		return tz;
	}

	if(Omega_M0 < 1.0 && Omega_Lambda0 == 0)
	{
		flag2 = 0;
		xx= step = 1.e-4;
		z1 = 1.0 + z;
		a = 1.0/z1;
		xxmax = 1.0e5;
		while(xx<xxmax && flag2 ==0)
		{
			atmp = 0.5*Omega_M0*(cosh(xx)-1.0)/(1.0-Omega_M0);
			if(atmp > a)
			{
				xxz = xx;
				flag2 = 1;
			}
			xx += step;
		}
		if(atmp < a)
		{
			printf("in hubble_time(z): increase your xxmax!\n Program exits!\n");
			exit(10);
		}
		tz = 0.5*Omega_M0*(sinh(xxz)-xxz)/pow(1.0-Omega_M0,1.5);
	//printf("f=2\n");
		return tz;
	}

	if(fabs(1.0-Omega_M0-Omega_Lambda0) < 1.e-4)
	{
		z1 = 1.0+z;
		omega1 = 1.0-Omega_M0;
		const0 = 2./(3.*sqrt(omega1));
		param0 = sqrt(omega1/Omega_M0)*pow(z1,-1.5);
		tz = const0*log(param0+sqrt(1.0+param0*param0));
	//printf("f=3\n");
		return tz;
	}

	printf("hubble_time(z) is not defined for this cosmology!\nProgram exits!\n");
	exit(10);
}


double lookbacktime(double z)
// computes lookbacktime in Gyrs at redshift z.
// lookbacktime is defined as
// lookbacktime = 0 at z=0
// lookbacktime = t0 at z=infty
{
	double lbt;

	lbt = (htime_z0 - hubble_time(z)) * xH0_recip;
	return lbt;
}

double z_hubble_time(double t1)
{
	double t, h0, a, z;
	t=t1*xhubble/10.;
	h0=1.02278 ; //H0=100km/s/Mpc
	a=pow(Omega_M0/Omega_Lambda0, 1.0/3)*pow(sinh(3./2*sqrt(Omega_Lambda0)*h0*t), 2./3);
	z=1./a-1.;
	return z;
}

double sigma_fit(double r)  
// copy Frank's program 
{
	double a1, a2, a3, a4, a5;
    double sigma_0, sigma;

	a1 = 71.59512;
    a2 = 1.23491;
//	a2=1.205; // this provides better fit
    a3 = -1.89001;
    a4 = 1.17110;
    a5 = -0.23848;

    sigma_0 = (1.0 + a2*pow(r,0.3) + a3*pow(r,0.4) +
                a4*pow(r,0.5) + a5*pow(r,0.6));
    sigma = a1*pow(sigma_0, -10.0);

    return sigma;
}


double var_Mo(double m) 
// copy Frank's program. [m]= h^-1 Msun 
{

        double rrr, r1, ri, r8, var_Mo;

        rrr = m/rho_aver0*3.0/4.0/PI;
        r1 = pow(rrr, 1.0/3.0);

        ri = r1 * Gamma_cosmo/0.25;
        r8 = 8.0 * Gamma_cosmo/0.25;

        var_Mo = sigma_fit(ri)*Sigma8/sigma_fit(r8);
//printf("%g %g %g\n",r1,m,var_Mo*var_Mo);
//printf("%g %g %g %g %g %g %g\n",rho_crit_0,Omega_M0,rrr,r1,ri,r8,var_Mo);
        return var_Mo;
}

double sig_Prada(double m)
// From Prada et al. 2012; a fitting formula for p(k) used for the Bolshoi simulation
{
	double y, sig;
	y = 1./(m/1e12);
	sig = 16.9 * pow(y, 0.41) / ( 1 + 1.102 * pow(y, 0.2) + 6.22 * pow(y, 0.333));
	return sig;
}


double transferWDM(double xk)
// the WDM transfer function T_WDM(k) = sqrt(P_WDM(k)/P_CDM(k))
// See Sommer-Larsen & Dolgov, 2001, ApJ, 551, 608
{
	double rf, fac, tr, tf;

	rf = 0.065 * pow(Omega_M0 -Omega_Baryon0, -1.0/3.0) * pow(MF_WDM/1.0e11, 1.0/3.0);

	fac = (xk*rf) + pow(xk*rf,2.0);

	tr = exp(-fac/2.0);

	return tr;
}


double power_spec(double xk)
// The CDM power spectrum with arbitrary normalization.
// The transfer function is taken from BBKS (Bardeen et el 1986) or
// from Efstathiou, Bond & White (1992) or from Eisenstein & Hu.
// The initial power-spectrum directly after inflation is a power-law
// with index `nspec'. For nspec=1 this yields the standard
// Harrison-Zel'dovich spectrum. The normalization is set by xk0
// which defines the wavelength at which the amplitude of the initial
// fluctuation spectrum is unity (this is arbitrary). The actual
// normalisation of the power-spectrum is set by sigma8.
{
        double xk0, tk, q, t1, t2, t3;
        double silk, stilde, fff, C1, C2, t11, t12, tb1, tb2;
        double T_c, T_b, T_WDM, ps;
// set normalization of initial power spectrum
        xk0 = 1./3000.0;
// initialize tk to check that computation was succesful
        tk = -1.0;
// the BBKS fitting function
        if(BBSK)
        {
                q = xk/Gamma_cosmo;
                t1 = log(1.0+2.34*q)/(2.34*q);
                t2 = 1.0 + (3.89*q) + pow(16.1*q, 2) + pow(5.46*q, 3) + pow(6.71*q, 4);
                tk = t1*pow(t2, -0.25);
        }

// the EBW fitting function
        if(EBW)
        {
                q = xk/Gamma_cosmo;
                tk = 1.0 + pow((6.4*q) + pow(3.0*q, 3.0/2) + pow(1.7*q, 2), 1.13);
                tk = pow(tk, -1.0/1.13);
        }

// the Eisenstein & Hu fitting function
        if(EISHU)
        {
                q = xk/(13.41*keq);

                silk = pow(xk/ksilk,1.4);
                silk = exp(-silk);
                stilde = sEH/pow(1.+pow(bnode/(xk*sEH), 3),0.3333);

                fff = 1.0/(1.0 + pow(xk*sEH/5.4, 4));
                C1 = 14.2 + 386.0/(1.0+69.0*pow(q,1.08));
                C2 = 14.2/alpha_c + 386.0/(1.0 + 69.9*pow(q,1.08));

                t11 = log(E + 1.8*beta_c*q);
                t12 = log(E + 1.8*q);

                t1 = t11/(t11 + C1*q*q);
                t2 = t11/(t11 + C2*q*q);
                t3 = t12/(t12 + C1*q*q);

                tb1 = t3/( 1.0 + pow(xk*sEH/5.2, 2));
                tb2 = (alpha_b/(1.0 + pow(beta_b/(xk*sEH), 3))) *silk;

                T_c = fff*t1 + (1.-fff)*t2;
                T_b = (tb1 + tb2) * sin(xk*stilde)/(xk*stilde);

                tk = (Omega_Baryon0/Omega_M0)*T_b + (Omega_M0-Omega_Baryon0)/Omega_M0 * T_c;
        }


// check that a transfer function has been defined
        if(tk < 0.0)
        {
                printf("transfer function has problem!\n");
                exit(10);
        }

// in the case of WDM, filter this

        if(MF_WDM == 0.0) T_WDM = 1.0;
        else T_WDM = transferWDM(xk);

        ps = tk*tk*pow(xk/xk0, nspec) * T_WDM*T_WDM;
        return ps;
}


double intvar(double xk, void *params)
{
        double x, wf, pk, f;
        double rf = *(double *)params;

        switch(ifilter)
        {
                case 1:
                        x = xk * rf;
                        wf = 3.0 * (sin(x) - x*cos(x)) / (x*x*x);
                        break;
                case 2:
                        x = xk * rf;
                        wf = exp(-0.5 * x *x);
                        break;
                case 3:
                        wf = 1.0;
                        break;
                default:
                        printf("You did not choose an appropriate window function!\n");
                        printf("exit from intvar!\n");
                        exit(2);
        }

        pk = power_spec(xk);
        f = pk * wf*wf * xk*xk;
//printf("%g %g\n", xk, pk);
        return f;
}

double var_numerical(double m)
{
        double rf, var;
        double ss, error;

        gsl_integration_workspace *w
                = gsl_integration_workspace_alloc(1000);

        gsl_function F;
        F.function = &intvar;
        F.params = &rf;

        switch(ifilter)
        {
                case 1:
                        rf = pow((3.*m)/(4.*PI*rho_aver0), 1./3);
                        gsl_integration_qagiu(&F, 0, 0, 1.e-7, 1000, w, &ss, &error);
                        break;
                case 2:
                        rf = 1./sqrt(2.*PI) * pow(m/rho_aver0, 1./3);
                        gsl_integration_qagiu(&F, 0, 0, 1.e-7, 1000, w, &ss, &error);
                        break;
                case 3:
                        rf = pow(m/(6.*PI*PI*rho_aver0), 1./3);
                        gsl_integration_qags(&F, 0, 1./rf, 0, 1.e-7, 1000, w, &ss, &error);
                        break;
                default:
                        printf("You did not choose an appropriate window function!\n");
                        exit(2);
        }
        var = Sigma8/sigma8_norm * sqrt(ss/(2.*PI*PI));
        return var;
}




