#define PI 3.14159265358979
#define E 2.7182818
#define YEAR_IN_SECOND 3.156e7
#define PC_IN_CM 3.0857e18
#define MSUN_IN_G 1.98892e33

extern double xH0;
extern double xH0_recip;
extern double htime_z0;
extern double H0;
extern double xhubble;
extern double Omega_M0;
extern double Omega_Lambda0;
extern double Omega_Baryon0_h2;
extern double Omega_Baryon0;
extern double MF_WDM;
extern double G;
extern double Sigma8;
extern double nspec;
extern int BBSK;
extern int EBW;
extern int EISHU;
extern int ifilter;
extern int ivar;
extern double Gamma_cosmo;
extern double rho_crit0;
extern double rho_aver0;


void init_cosmo(void);
double xH(double z); 
double Omega_m(double z); 
double rho_crit(double z); // critical density of the universe
double lin_growth(double z); 
double delta_c(double z);
double Delta_vir(double z);
double hubble_time(double z);
double lookbacktime(double z);
double z_hubble_time(double t1);
double sigma_fit(double r);  
double var_Mo(double m); 
double sig_Prada(double m);
double transferWDM(double xk);
double power_spec(double xk);
double intvar(double xk, void *params);
double var_numerical(double m);

