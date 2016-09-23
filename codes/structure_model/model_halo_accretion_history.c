/*
 * model_halo_accretion_history.c
 *
 *  Created on: Jun 14, 2012
 *      Author: luyu
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "variables.h"
#include "proto.h"
#include "cosmo.h"

double mah_wechsler(double a)
{
	double s=2,a0=1,alpha, m0 = pow(10., 11);
	double ac,m;

	alpha=1.65;
	ac = a0*alpha/s;
//ac=0.75; // for 1.e14msun halo
//ac=0.6; // for 1.e13msun halo
//ac=0.53; // for 1.e12msun halo
ac=0.45; // for 1.e11msun halo
//ac = 0.50;
	m = m0 * exp(-ac/a0*s*(a0/a - 1));
	//m = exp(-s*0.15873/a );
	return m;
}

double mah_mcbride(double z)
{
	double beta, gamma, m0, m;
	double lgm0;

	lgm0 = Mass_bin;
	beta = 0.8;
	//printf("lgm0=%g, beta=%g\n", lgm0, beta);
	/*
	if (Mass_bin == 12) {beta = 0.6; lgm0 = 12;}
	if (Mass_bin == 11) {beta = 0.8; lgm0 = 11;}
	if (Mass_bin == 12.1) {beta = 0.6; lgm0 = 12.1;}
	if (Mass_bin == 11.1) {beta = 0.8; lgm0 = 11.1;}
	*/
	if (Mass_bin > 11) beta = (0.6-0.8) / (12.-11) * (lgm0-11.) + 0.8;

	gamma = 0.9;
	m0 = 1*pow(10., lgm0);

	m = m0*pow(1+z, beta) / exp(gamma*z) * exp(gamma * Redshift_end)/pow(1+Redshift_end, beta);
	return m;
}

double mah(double z)
{
	double m;
	double a = 1./(1+z);
	if (Mah_simu)
	{
		m = mah_simu(z);
		/*
		printf("m=%g\n",m);
		if (m <=0.0) m = mah_simu(0.0) * mah_wechsler(a)/mah_wechsler(1);
		*/
	}
	else
	{
	 //m = mah_wechsler(a);
	 m = mah_mcbride(z);
	}

	return m;
}

double mhdot_ma(double z, double mh)
{
	double a, a1, a2, z1, z2, t1, t2, m1, m2;
    double y;

    a = 1./(1+z);
    a1 = a-0.001;
    a2 = a+0.001;
    z1 = 1./a1 - 1;
    z2 = 1./a2 - 1;
    t1 = hubble_time(z1) * xH0_recip;
    t2 = hubble_time(z2) * xH0_recip;

    m1 = mah(z1);
    m2 = mah(z2);

    y = (m2-m1)/(t2-t1);

    return y;
}

double mhdot_dekel(double z, double mh)
{
	double mhd;
	//mhd = 0.47*mh*pow(mh/1.e12,0.15)*pow((1+z)/3,2.25);
	mhd = 34.*1e9 * pow(mh/1e12, 1.14)*pow((1.+z), 2.4);
	return mhd;
}

double halo_mass_accretion_rate(double mvir, double z)
{
	//return mhdot_dekel(z, mvir);
	//printf("hma_rate: %g %g\n", z, mvir);
	return mhdot_ma(z, mvir);
}

double mass_profile_nfw(double x, double c)
{
	double m0, mr;

	m0 = 1.0/(log(1.+c) - c/(1.+c));
	mr = m0 *(log(1.+x*c) - x*c/(1.+x*c));
	return mr;
}

void halo_mass_profile(struct galaxy *gal)
{
	int i;
	double r;
//printf("halo_mass_profile: c=%g\n", gal->ConcenHalo);
	for(i=0; i<gal->nbin; i++)
	{
		r = gal->RadiusInner[i] / gal->RadiusHalo;
		if (r <=1) gal->MassProfHalo[i] = gal->MassHalo * mass_profile_nfw(r, gal->ConcenHalo);
		else gal->MassProfHalo[i] = gal->MassHalo;
	}
}

double fun_cmin(double x)
{
	double c0 = 3.681, c1 = 5.033, alpha = 6.948, x0 = 0.424;
	double r;

	r = c0 + (c1-c0)*(1./M_PI * atan(alpha * (x-x0))+0.5);
	return r;
}

double fun_smin(double x)
{
	double s0 = 1.047, s1 = 1.646, beta = 7.386, x1 = 0.526;
	double r;

	r = s0 + (s1 - s0) * (1./M_PI * atan(beta*(x-x1))+0.5);
	return r;
}

double halo_concentration_prada(double z, double m)
// model for c_200
{
	double ap = 2.881, bp = 1.257, cp = 1.022, dp = 0.060;
	double x, mh, d, sig, b0, b1, sig_p, cterm, c;

	x = pow(Omega_Lambda0/Omega_M0, 1./3);
	mh = m * xhubble; // the rest of the calculation is in unit of Msun/h;
	d = lin_growth(z);
	sig = d * sig_Prada(m);
	b0 = fun_cmin(x)/fun_cmin(1.393);
	b1 = fun_smin(x)/fun_smin(1.393);

	sig_p = b1 * sig;
	cterm = ap * (pow(sig_p/bp, cp)+1) * exp(dp / sig_p / sig_p);

	c = b0 * cterm;
	return c;
}

double halo_concentration(double z, struct galaxy *gal)
{
	double c;

	//c= (1.+2.3)/(1+z) * 3.5;
	c = 1.35 * halo_concentration_prada(z, gal->MassHalo); //convert into c_vir from c_200
	return c;
}

double halo_vmax(double vvir, double c)
{
	double v2, vmax2, a0, vmax;

	v2 = vvir*vvir;
	a0 = 1.0/(log(1.+c) - c/(1.+c));
	vmax2 = 0.216 * c / a0 * v2;
	vmax = sqrt(vmax2);

	return vmax;
}

int mah_nz, mah_nhalo;
double *mah_z, *mah_a, *mah_lgm, **mah_array;

int read_simu_mah(void)
{
	int i, j, nhalo, nz;
	char filename[200], sbuf[200];
	FILE *fd;

	//sprintf(filename,"./mah_180_m12.0.txt");

    // NATHAN'S MAH DIRECTORIES
     if (Redshift_end == 0.0) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_less2/mah_180_m10_12_selected.txt");
     if (Redshift_end == 0.5) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_less2/mah_094_m10_12_selected.txt");
     if (Redshift_end == 1.0) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_less2/mah_066_m10_12_selected.txt");
     if (Redshift_end == 2.0) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_less2/mah_039_m10_12_selected.txt");
     if (Redshift_end == 3.0) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_less2/mah_027_m10_12_selected.txt");
     if (Redshift_end == 4.0) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_less2/mah_019_m10_12_selected.txt");
     if (Redshift_end == 5.0) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_less2/mah_013_m10_12_selected.txt");
     if (Redshift_end == 6.0) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_less2/mah_009_m10_12_selected.txt");
     
     if (Redshift_end == 0.0) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_chin/mah_099_m10_12_selected.txt");  // Chinchilla simulation
     if (Redshift_end == 0.5) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_chin/mah_084_m10_12_selected.txt");  // Chinchilla simulation
     if (Redshift_end == 1.0) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_chin/mah_073_m10_12_selected.txt");  // Chinchilla simulation
     if (Redshift_end == 2.0) sprintf(filename, "/Users/Nathan/Documents/Carnegie/Summer_2015/Github/summer_project2015/mah/mah_chin/mah_058_m10_12_selected.txt");  // Chinchilla simulation
    

    /* YU'S MAH DIRECTORIES
	if (Redshift_end == 0.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_180_m10_12_selected.txt");
	if (Redshift_end == 0.5) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_094_m10_12_selected.txt");
	if (Redshift_end == 1.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_066_m10_12_selected.txt");
	if (Redshift_end == 2.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_039_m10_12_selected.txt");
	if (Redshift_end == 3.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_027_m10_12_selected.txt");
	if (Redshift_end == 4.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_019_m10_12_selected.txt");
	if (Redshift_end == 5.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_013_m10_12_selected.txt");
	if (Redshift_end == 6.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_009_m10_12_selected.txt");

	if (Redshift_end == 0.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_chin/mah_099_m10_12_selected.txt");  // Chinchilla simulation
	if (Redshift_end == 0.5) sprintf(filename, "/Users/luyu/project/disc/results/mah_chin/mah_084_m10_12_selected.txt");  // Chinchilla simulation
	if (Redshift_end == 1.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_chin/mah_073_m10_12_selected.txt");  // Chinchilla simulation
	if (Redshift_end == 2.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_chin/mah_058_m10_12_selected.txt");  // Chinchilla simulation
    */

/*
  // mass function weighting
	if (Redshift_end == 0.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_180_m10_14_selected_0.1.txt");
	if (Redshift_end == 0.5) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_094_m10_14_selected_0.1.txt");
	if (Redshift_end == 1.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_066_m10_14_selected_0.1.txt");
	if (Redshift_end == 2.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_039_m10_14_selected_0.1.txt");
	if (Redshift_end == 3.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_027_m10_14_selected_0.1.txt");
	if (Redshift_end == 4.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_019_m10_14_selected_0.1.txt");
	if (Redshift_end == 5.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_013_m10_14_selected_0.1.txt");
	if (Redshift_end == 6.0) sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_009_m10_14_selected_0.1.txt");
*/
	//sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_selected.txt");

	if(Mass_bin == 11.0)
	{
			//sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_180_m11.0.txt");
		sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_180_m11.0_selected.txt");
	}
	if(Mass_bin == 12.0)
	{
		sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_180_m12.0.txt");
	}
	if(Mass_bin > 11.0 && Mass_bin <11.5)
	{
		//sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_180_m11.1.txt");
		sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_180_m11.1_selected.txt");
	}
	if(Mass_bin > 12.0 && Mass_bin <12.5)
	{
		printf("yes!\n");
		sprintf(filename, "/Users/luyu/project/disc/results/mah_less2/mah_180_m12.1.txt");
	}
	if(!(  fd=fopen(filename,"r")))
	{
		printf("I cannot open the file '%s'. %g %g\n",filename, Redshift_end, Mass_bin);
		exit(1);
	}
	//fscanf(fd, "%s", sbuf);
	fgets(sbuf, 200, fd);
	fscanf(fd, "%d %d", &nhalo, &nz);
	printf("read_simu_mah: %s %d %d\n", sbuf, nhalo, nz);
	mah_nz = nz;
	mah_nhalo = nhalo;
	mah_z=(double *)malloc(nz * sizeof(double));
	mah_lgm = (double *)malloc(nz * sizeof(double));
	for(j=0; j<nz; j++)
		fscanf(fd, "%lf", &mah_z[j]);
	mah_array=(double **)malloc(nhalo*sizeof(double *));
	for (i = 0; i < nhalo; i++)
	{
	    mah_array[i] = (double *) malloc( nz * sizeof(double));
	    for(j=0; j < nz; j++)
	    {
	    	fscanf(fd, "%lf", &mah_array[i][j]);
	    }
	    //printf("read_simu_mah: %d\n", i);
	}
	fclose(fd);
	//exit(0);
	return nhalo;
}

void select_simu_mah(int iseed)
{
	int i;
	int ihalo;

	if(iseed < 0) ihalo = (int) random_number_gstructure(0)*(mah_nhalo-1);
	else ihalo = iseed;
	if( ihalo >= mah_nhalo)
	{
		printf("not that many halos in the MAH file! using the first halo instead!");
		ihalo = 0;
	}

	for(i=0;i<mah_nz;i++)
	{
		mah_lgm[i] = dmax(log10(mah_array[ihalo][i]) + 10., 6);
		if(i>1) mah_lgm[i] = dmax(mah_lgm[i], mah_lgm[i-1]);
		//printf("select_simu_mah: %g %g\n", mah_z[i], mah_lgm[i]);
	}
	//exit(0);
}

void free_simu_mah(void)
{
	int i;
	free(mah_z);
	free(mah_a);
	free(mah_lgm);
	for(i = 0; i<mah_nhalo; i++)
		free(mah_array[i]);
	free(mah_array);
}

double mah_simu(double z)
{
	double lgm, m;

	if(z == 0)
		lgm = mah_lgm[mah_nz-1];
	else
		lgm = interpolate_bipoint(mah_z, mah_lgm, mah_nz, z, -1);
	m = pow(10, lgm);

	//printf("mah_simu: %g %g\n", z, m);
	//m = mah_mcbride(z);
	return m;
}
