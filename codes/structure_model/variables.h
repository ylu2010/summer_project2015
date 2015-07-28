
#define N_RADIUS_BIN 200
#define N_TIME_BIN 500

#define BOLTZMANN   1.3806e-16
#define PROTONMASS  1.6726e-24
#define ELECTRONVOLT 1.60217653e-12  // erg


extern double UnitMass_in_g;
extern double UnitLength_in_cm;
extern double UnitTime_in_Megayears;
extern double UnitTime_in_Second;
extern double UnitVelocity_in_cm_per_s;
extern double BaryonFrac;
extern double SolarMetallicity;
extern double MinimumMetallicityRelativeToSolar;

extern int Metal_gas_evolu;
extern int Mah_simu;
extern int Do_preheating;
extern int Do_reinfall;
extern float Mass_bin;
extern double Redshift;
extern double Redshift_end;
extern double Bin_size_time;

struct parameter
{
	double StarFormationEfficiency;
	double DiskRadiusFactor;
	double StarFormationCriticalSurfaceDensity;
	double SNLoadingFactor;
	double SNLoadingFactorIndex;
	double BaryonAccretionFraction;
	double StellarMassLossFraction;
	double Reionization_z0;
	double Reionization_zr;
	double PreheatEntropy;
	double PreheatEntropySlope;
	double EntropyRatio;
	double EntropyProfileIndex;
	double Yield;
	double ZFractionYieldToEject;
	double ZFractionYieldToHot;
	double MassFractionEjectToHot;
	double GalaxyHeatingEfficiency;
};

struct galaxy
{
	int nbin;
	int ntbin; // time array
    double z;
    double MassBin; //new //to better separate halos in data tables
    double MassHalo;
    double MassHot;
    double MassCloud;
    double MassCold;
    double MassStar;
    double MassEject;
    double MassColdMolecular;
    double MassColdAtomic;
    double MassColdIonized;
    double MetalHot;
    double MetalCold;
    double MetalStar;
    double MetalEject;
    double RadiusHalo;
    double ConcenHalo;
    double SpinHalo;
    double SpinCooling;
    double RadiusCooling;
    double RadiusDisc;
    double RadiusHalfStar;
    double RadiusHalfCold;
    double RateHaloAccretion;
    double VelocityVirial;
    double VelocityMax;
    double RateCooling;
    double RateStarFormation;
    double RateOutflow;
    double TemperatureVirial;
    double EntropyVirial;
    double Entropy;
    double TimeCooling;
    double RadiusMostInner;
    double RadiusMostOuter;
    double RadiusInner[N_RADIUS_BIN];
    double RadiusOuter[N_RADIUS_BIN];
    double MassProfHalo[N_RADIUS_BIN];
    double MassProfDM[N_RADIUS_BIN];
    double MassProfDMContracted[N_RADIUS_BIN];
    double MassProfHot[N_RADIUS_BIN];
    double MassProfStar[N_RADIUS_BIN];
    double MassProfCold[N_RADIUS_BIN];
    double MassProfNeutral[N_RADIUS_BIN];
    double DensityProfHot[N_RADIUS_BIN];
    double TemperatureProfHot[N_RADIUS_BIN];
    double CoolingRate[N_RADIUS_BIN];
    double CoolingTime[N_RADIUS_BIN];
  //double MassMetalHot[N_RADIUS_BIN]; //not tracked; not in print_galaxy
    double MassMetalCold[N_RADIUS_BIN];
    double MassMetalStar[N_RADIUS_BIN];
    double MetallicityCold[N_RADIUS_BIN]; //new
    double MetallicityStar[N_RADIUS_BIN]; //new
    double SDensityMetalCold[N_RADIUS_BIN]; //new
    double SDensityMetalStar[N_RADIUS_BIN]; //new
    double SDensityCold[N_RADIUS_BIN];
    double SDensityColdMolecular[N_RADIUS_BIN];
    double SDensityColdAtomic[N_RADIUS_BIN];
    double SDensityStar[N_RADIUS_BIN];
    double SDensitySFR[N_RADIUS_BIN];
	double SDensityOFR[N_RADIUS_BIN];
    double SDensityCAR[N_RADIUS_BIN];//Cold gas accretion rate
    double TimeArray[N_TIME_BIN];
    double StarFormationHistory[N_TIME_BIN]; // stellar mass formed in a time interval
    double SDensitySFH[N_RADIUS_BIN][N_TIME_BIN];
};

extern struct parameter Par;
extern FILE *fp_hist;
extern FILE *fp_disc;
extern FILE *fp_pred;
extern FILE *fp_list;
extern FILE *fp_snap;

extern int Write_pred_file;
extern int Write_hist_file;
extern int Write_prof_file;
extern int Write_snap_file;

struct interval
{
		double min, max;
};
