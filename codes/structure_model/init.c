/*
 * init.c
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

double Bin_size_time;

void init(struct galaxy *gal)
{
	int i, j, nbin=N_RADIUS_BIN, ntbin=N_TIME_BIN;
	double lgrmin = -4.0, dlgr;
	dlgr = (0-lgrmin)/nbin;

	gal->nbin = nbin;
	gal->ntbin = ntbin;
    gal->z = 100;
    gal->MassHalo = 0.0;
    gal->MassHot = 0.0;
    gal->MassCloud=0.0;
    gal->MassCold = 0.0;
    gal->MassStar = 0.0;
    gal->MassColdAtomic = 0.0;
    gal->MassColdMolecular = 0.0;
    gal->MassColdIonized = 0.0;
    gal->MassEject = 0;
    gal->MetalHot = 0.0;
    gal->MetalCold = 0.00;
    gal->MetalStar = 0.0;
    gal->MetalEject = 0.0;
    gal->ConcenHalo = 10.0;
    gal->RadiusHalo = 0.0;
    gal->RadiusDisc = 0.00;
    gal->RadiusHalfCold = 0.0;
    gal->RadiusHalfStar = 0.0;
    gal->RateStarFormation = 0.0;
    gal->VelocityVirial = 0.0;
    gal->SpinHalo = 0.035;
    gal->SpinCooling = 1;
    gal->TimeCooling = 0.0;
    /*
    gal->RadiusInner = malloc(gal->nbin * sizeof(double));
    gal->RadiusOuter = malloc(gal->nbin * sizeof(double));
    gal->MassProfHalo = malloc(gal->nbin * sinzeof(double));
    gal->SDensityCold = malloc(gal->nbin * sizeof(double));
    gal->SDensityStar = malloc(gal->nbin * sizeof(double));
	*/
    gal->RadiusMostOuter = 0.270; // 200kpc
    gal->RadiusMostInner = pow(10., lgrmin) * gal->RadiusMostOuter;
    for (i=0; i<nbin; i++)
    {
    	gal->RadiusInner[i] = pow(10., lgrmin+i*dlgr) * gal->RadiusMostOuter;
    	gal->RadiusOuter[i] = pow(10., lgrmin+(i+1)*dlgr) * gal->RadiusMostOuter;
    	gal->SDensityCold[i] = 0.0;
    	gal->SDensityStar[i] = 0.0;
    	gal->SDensityColdAtomic[i] = 0.0;
    	gal->SDensityColdMolecular[i] = 0.0;
	gal->MassMetalCold[i] = 0.0;
	gal->MassMetalStar[i] = 0.0;
    	gal->CoolingTime[i]=0.0;
    	for(j=0; j<ntbin; j++)
    	{
    		gal->SDensitySFH[i][j] = 0.0;
    	}
    }
    Bin_size_time = hubble_time(0.0)/ntbin * xH0_recip;
    for(j=0; j<ntbin; j++)
    {
    	gal->TimeArray[j] = Bin_size_time * (j+0.5);
    	gal->StarFormationHistory[j] = 0.0;
    	//printf("init: j=%d t=%g\n", j, gal->TimeArray[j]);
    }

}
