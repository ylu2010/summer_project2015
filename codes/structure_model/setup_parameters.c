/*
 * setup_parameters.c
 *
 *  Created on: September 21, 2016
 *      Author: luyu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "variables.h"
#include "proto.h"
//#include "ihs.hpp"

#define DOUBLE 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

double Reionization_z0;
double Reionization_zr;
double BaryonAccretionFraction;
double StellarMassLossFraction;
double Yield;
double MassFractionEjectToHot;
double ReincorporationTimeScale;
double PreventionMassScale;
double PreventionMassIndex;
double PreventionRedshift;
double PreheatEntropySlope;
double PreheatEntropy;
double EntropyProfileIndex;
double DiskRadiusFactor;
double ZFractionYieldToEject;
double ZFractionYieldToHot;
double GalaxyHeatingEfficiency;
double SNLoadingFactor;
double SNLoadingFactorIndex;
double ReincorporationTimeScale;
double StarFormationCriticalSurfaceDensity; //Msun/pc^2 no use for the Krumholz model
double StarFormationEfficiency; // no use for Krumholz model
double EntropyRatio;

char Cooling_function_dir[400];
char Mah_file_name[400];
char OutputDir[400];

void read_parameter_file(char *fname)
{
	int i, j, nt=0;
	char buf[400], buf1[400], buf2[400], buf3[400];
	int id[MAXTAGS];
	void *addr[MAXTAGS];
	char tag[MAXTAGS][50];
	int errorFlag = 0;
	FILE *fd;

    printf("\nreading parameter file: %s!\n\n", fname);
    
	strcpy(tag[nt], "OutputDir");
	addr[nt] = OutputDir;
	id[nt++] = STRING;

	strcpy(tag[nt], "Cooling_function_dir");
	addr[nt] = Cooling_function_dir;
	id[nt++] = STRING;

    strcpy(tag[nt], "Mah_simu");
    addr[nt] = &Mah_simu;
    id[nt++] = INT;
    
    strcpy(tag[nt], "Mah_file_name");
    addr[nt] = Mah_file_name;
    id[nt++] = STRING;

    strcpy(tag[nt], "Mass_bin");
    addr[nt] = &Mass_bin;
    id[nt++] = DOUBLE;
    
    strcpy(tag[nt], "BaryonFrac");
    addr[nt] = &BaryonFrac;
    id[nt++] = DOUBLE;
    
    strcpy(tag[nt], "Metal_gas_evolu");
    addr[nt] = &Metal_gas_evolu;
    id[nt++] = INT;
    
    strcpy(tag[nt], "SolarMetallicity");
    addr[nt] = &SolarMetallicity;
    id[nt++] = DOUBLE;
    
    strcpy(tag[nt], "MinimumMetallicityRelativeToSolar");
    addr[nt] = &MinimumMetallicityRelativeToSolar;
    id[nt++] = DOUBLE;
    
    strcpy(tag[nt], "Do_preheating");
    addr[nt] = &Do_preheating;
    id[nt++] = INT;
    
    strcpy(tag[nt], "Do_reinfall");
    addr[nt] = &Do_reinfall;
    id[nt++] = INT;

	strcpy(tag[nt], "Star_formation_model");
	addr[nt] = &Star_formation_model;
	id[nt++] = INT;
    
    strcpy(tag[nt], "Redshift_end");
    addr[nt] = &Redshift_end;
    id[nt++] = DOUBLE;
    
	strcpy(tag[nt], "Write_pred_file");
    addr[nt] = &Write_pred_file;
    id[nt++] = INT;

	strcpy(tag[nt], "Write_pred_saparately");
    addr[nt] = &Write_pred_saparately;
    id[nt++] = INT;

	strcpy(tag[nt], "Write_hist_file");
    addr[nt] = &Write_hist_file;
    id[nt++] = INT;

	strcpy(tag[nt], "Write_prof_file");
    addr[nt] = &Write_prof_file;
    id[nt++] = INT;

	strcpy(tag[nt], "Write_snap_file");
    addr[nt] = &Write_snap_file;
    id[nt++] = INT;

    strcpy(tag[nt], "Reionization_z0");
    addr[nt] = &Reionization_z0;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "Reionization_zr");
    addr[nt] = &Reionization_zr;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "BaryonAccretionFraction");
    addr[nt] = &BaryonAccretionFraction;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "StellarMassLossFraction");
    addr[nt] = &StellarMassLossFraction;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "Yield");
    addr[nt] = &Yield;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "MassFractionEjectToHot");
    addr[nt] = &MassFractionEjectToHot;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "ReincorporationTimeScale");
    addr[nt] = &ReincorporationTimeScale;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "PreventionMassIndex");
    addr[nt] = &PreventionMassIndex;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "PreheatEntropySlope");
    addr[nt] = &PreheatEntropySlope;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "PreheatEntropy");
    addr[nt] = &PreheatEntropy;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "EntropyProfileIndex");
    addr[nt] = &EntropyProfileIndex;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "DiskRadiusFactor");
    addr[nt] = &DiskRadiusFactor;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "ZFractionYieldToEject");
    addr[nt] = &ZFractionYieldToEject;
    id[nt++] = DOUBLE;

    strcpy(tag[nt], "ZFractionYieldToHot");
    addr[nt] = &ZFractionYieldToHot;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "GalaxyHeatingEfficiency");
    addr[nt] = &GalaxyHeatingEfficiency;
    id[nt++] = DOUBLE;
    
	strcpy(tag[nt], "SNLoadingFactor");
    addr[nt] = &SNLoadingFactor;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "SNLoadingFactorIndex");
    addr[nt] = &SNLoadingFactorIndex;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "StarFormationCriticalSurfaceDensity");
    addr[nt] = &StarFormationCriticalSurfaceDensity;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "StarFormationEfficiency");
    addr[nt] = &StarFormationEfficiency;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "EntropyRatio");
    addr[nt] = &EntropyRatio;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "PreventionMassScale");
    addr[nt] = &PreventionMassScale;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "PreventionMassIndex");
    addr[nt] = &PreventionMassIndex;
    id[nt++] = DOUBLE;

	strcpy(tag[nt], "PreventionRedshift");
    addr[nt] = &PreventionRedshift;
    id[nt++] = DOUBLE;


    if((fd = fopen(fname, "r")))
    {
        while(!feof(fd))
        {
            *buf = 0;
            fgets(buf, 200, fd);
            if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
                continue;
            
            if(buf1[0] == '%')
                continue;
            
            for(i = 0, j = -1; i < nt; i++)
                if(strcmp(buf1, tag[i]) == 0)
                {
                    j = i;
                    tag[i][0] = 0;
                    break;
                }
            
            if(j >= 0)
            {
                printf("%35s\t%10s\n", buf1, buf2);
                
                switch (id[j])
                {
                    case DOUBLE:
                        *((double *) addr[j]) = atof(buf2);
                        break;
                    case STRING:
                        strcpy(addr[j], buf2);
                        break;
                    case INT:
                        *((int *) addr[j]) = atoi(buf2);
                        break;
                }
            }
            else
            {
                printf("Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname, buf1);
                errorFlag = 1;
            }
        }
        fclose(fd);
        
        i = strlen(OutputDir);
        if(i > 0)
            if(OutputDir[i - 1] != '/')
                strcat(OutputDir, "/");
    }
    else
    {
        printf("Parameter file %s not found.\n", fname);
        errorFlag = 1;
    }

    for(i = 0; i < nt; i++)
    {
        if(*tag[i])
        {
            printf("Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
            errorFlag = 1;
        }
    }
    
    if(errorFlag)
        exit(1);

	printf(" Mass_bin = %f\n", Mass_bin);
	printf(" BaryonFrac = %f\n", BaryonFrac);
}
	

void init_parameters(void)
{
	/* print out run constants */
	printf("You are running the model with:\n");
	printf(" BaryonFrac = %f\n", BaryonFrac);
	printf(" SolarMetallicity = %f\n", SolarMetallicity);
	printf(" MinimumMetallicityRelativeToSolar = %f\n", MinimumMetallicityRelativeToSolar);
	printf(" Metal_gas_evolu : %d\n", Metal_gas_evolu);
	printf(" Mah_simu : %d\n", Mah_simu);
	printf(" Do_preheating : %d\n", Do_preheating);
	printf(" Star_formation_model : %d\n", Star_formation_model);
	printf(" Do_reinfall = %d\n", Do_reinfall);
	printf(" Mass_bin = %f\n", Mass_bin);
	printf(" Write_pred_file : %d\n", Write_pred_file);
	printf(" Write_pred_saparately : %d\n", Write_pred_saparately);
	printf(" Write_hist_file : %d\n", Write_hist_file);
	printf(" Write_prof_file : %d\n", Write_prof_file);
	printf(" Write_snap_file : %d\n", Write_snap_file);
	printf(" Redshift_end = %f\n", Redshift_end);

	/* setup run parameters */
    Par.Reionization_z0 = Reionization_z0;
    Par.Reionization_zr = Reionization_zr;

    Par.BaryonAccretionFraction = BaryonAccretionFraction;
    Par.StellarMassLossFraction = StellarMassLossFraction; 
    Par.Yield = Yield;

    Par.MassFractionEjectToHot = MassFractionEjectToHot;
    Par.ReincorporationTimeScale = ReincorporationTimeScale;

	Par.PreventionMassScale = PreventionMassScale;
    Par.PreventionMassIndex = PreventionMassIndex;
	Par.PreventionRedshift = PreventionRedshift;

    Par.PreheatEntropySlope = PreheatEntropySlope;
    Par.PreheatEntropy = PreheatEntropy;
    Par.EntropyProfileIndex = EntropyProfileIndex;
    Par.DiskRadiusFactor = DiskRadiusFactor;
    Par.ZFractionYieldToEject = ZFractionYieldToEject;
    Par.ZFractionYieldToHot = ZFractionYieldToHot; 
    Par.GalaxyHeatingEfficiency = GalaxyHeatingEfficiency;
    Par.SNLoadingFactor = SNLoadingFactor;
    Par.SNLoadingFactorIndex = SNLoadingFactorIndex;
    Par.ReincorporationTimeScale = ReincorporationTimeScale;

    Par.StarFormationCriticalSurfaceDensity = StarFormationCriticalSurfaceDensity; //Msun/pc^2 no use for the Krumholz model
    Par.StarFormationEfficiency = StarFormationEfficiency; // no use for Krumholz model
    Par.EntropyRatio = EntropyRatio;

}
