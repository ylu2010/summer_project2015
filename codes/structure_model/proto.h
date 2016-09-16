/*
 * proto.h
 *
 *  Created on: Jun 14, 2012
 *      Author: luyu
 */

#ifndef PROTO_H_
#define PROTO_H_
#endif /* PROTO_H_ */

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
extern "C" {
#endif

int set_varying_parameters(double *params, int nparams);
int set_predictions(double *p, int np, struct galaxy *gal, int ihalo);

// run_galaxy.c
int setup_run(void);
int run_galaxy(double *params, int nparams, double *preds, int npreds, int mode, int ip);
int finalize_run(void);


void init(struct galaxy *gal);
//void resize_radius_bins(struct galaxy *gal);

void cold_gas_accretion(struct galaxy *gal, double t, double dt);
void cold_gas_accretion_surface(struct galaxy *gal, double t, double dt);

void star_formation_global(struct galaxy *gal, double dt);

double molecular_fraction(double den, double z0);
void star_formation_surface(struct galaxy *gal, double t, double dt);
void star_formation_surface_molecule(struct galaxy *gal, double t, double dt);
void star_formation_surface_molecule_with_guo2011_feedback(struct galaxy *gal, double t, double dt);
void ejected_gas_reincorporation(struct galaxy *gal, double t, double dt);
//void star_formation_surface_limit(struct galaxy *gal, double dt);
double outflow_massloading_factor(struct galaxy *gal);

void star_formation_rate_molecule(struct galaxy *gal, double dt);
void disc_mass_composition(struct galaxy *gal);

void evolve_galaxy(struct galaxy *gal, int write_in_file);

// model_halo_accetion.c
double halo_mass_accretion_rate(double mvir, double z);
double mass_profile_nfw(double x, double c);
void halo_mass_profile(struct galaxy *gal);
double halo_concentration_prada(double z, double m);
double halo_concentration(double z, struct galaxy *gal);
double halo_vmax(double vvir, double c);
double mah_wechsler(double a);
double mah_mcbride(double z);
double mah(double z);
int read_simu_mah(void);
void select_simu_mah(int iseed);
void free_simu_mah(void);
double mah_simu(double z);

// model_hotgas_accretion.c
double do_reionization(double mvir, double Zcurr);
double do_preheating(struct galaxy *gal, double zcurr);
//double hot_gas_accretion_rate(struct galaxy *gal, double z);
double hot_gas_accretion_fraction(struct galaxy *gal, double z);
double halo_cold_gas_accretion_fraction(struct galaxy *gal, double z);

void hot_gas_profile(struct galaxy *gal);
void hot_gas_profile_isenthermal(struct galaxy *gal);
void hot_gas_profile_power_law_entropy ( struct galaxy *gal );

double disk_radius(struct galaxy *gal);

void init_file(void);
void close_file(void);

void sig_metal_calc(struct galaxy *gal);
void halo_adjust(struct galaxy *gal, double z, double mh);
void print_galaxy(struct galaxy *gal);
void print_snapshot(struct galaxy *gal, double z,double thubble, double dt, double mh);

double dmax(double x, double y);
double dmin(double x, double y);
double interpolate_bipoint(double *x, double *y, int n, double x0, int flag);
double second(void);
void init_random_number_gstructure(void);
void free_random_number_gstructure(void);
double random_number_gstructure(int flag);

// cooling.c
void read_cooling_function(void);
void free_cooling_table(void);
double cooling_rate(struct galaxy *gal, double t, double dt, int model);
double cooling_rate_shell(struct galaxy *gal, double hubble_time, double dt, int model);

// model_adiabatic_contraction.c
void adiabatic_contraction(struct galaxy *gal);

#ifdef __cplusplus /* If this is a C++ compiler, use C linkage */
}
#endif
