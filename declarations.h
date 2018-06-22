#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_spline.h>

#define Z_0       0.3  // peak of n(z) is at 2*Z_0

#define DO_NONLIN 1

#define DELTA_C   1.686 /* VERY IMPORTANT - 1.686, not 1.69 */
#define DLNM    0.05

#define SMALL  1e-6
#define VERYSMALL  1e-15

#define M_SUN 1.989e33                             /* grams       */

#define MPC_TO_CM  (3.0856e24)
#define CRIT_DENS   (1.8791e-29*hhubble*hhubble*pow(MPC_TO_CM, 3)) /* in grams Mpc^{-3} */

#define K_START_INTERP 2.0e-5 // make sure this is covered by output of CAMB
#define K_END_INTERP   1.0e2  // making sure this is covered by output of CAMB

#define LOG10_KMIN_CIJ_INT -3.0 // min range of C_ij integral in  k, CHECKED -3 IS GOOD FOR ALL z (Oct 2016)
#define LOG10_KMAX_CIJ_INT  0.5 // max range ; FOUND NEEDS TO BE AT LEAST =0.5 BECAUSE OF MAD OSC OF jl'(x) UP TO HIGH X (Oct 2016)

#define Z_EQ  3000.0  // doesnt matter exactly what it is, as long as its high-z 
/***********************************************/
/** l_range of completeness - 20 is fiducial  **/
/***********************************************/
#define L_MAX_CALIB 10

/*******************************************************************/
/** histogram spacing for bias/err for individual (l1, m1)        **/
/** large value needed for large, super-duper safe, range of hist **/
/**        right now, range from 1.0e-1 to 1.0e4                  **/
/*******************************************************************/
#define BINS_HIST     1000000
#define HIST_MAXVALUE 500.0

/*********************************************************************/
/**                l_range of data **/
/** make sure you got enough of these at l<l_nonlin for large L_MAX **/
/** notice you are doing every ell below CONSECUTIVE_ELL_LIMIT      **/
/*********************************************************************/
#define CONSECUTIVE_ELL_LIMIT 20  // do all ell consecutively from ell=2 up to here; 20 fiducial; 1000 to make S/N plot/
#define L_MAX    1000.0

#define L_BINS_ABOVE_CONSEC 10 // this is where you pick how many 'big bins'; 10 is fiducial; 0 to make S/N plot

// this is now in CL_LSS.c, and there's an option to increase it by 1 if starting with dipole
//#define L_BINS   (L_BINS_ABOVE_CONSEC + CONSECUTIVE_ELL_LIMIT-1)

#define DL    ((L_MAX-CONSECUTIVE_ELL_LIMIT) / (L_BINS_ABOVE_CONSEC-0.5) )

#define MIN_ELL_USE_APPROX 50 // min ell at which you start using the approx expression

/*****************************************************************************/
/** where you stop using FULL and start using LIMBER power spectrum         **/
/** 30.5 fails with kr<100*ell (when running SN_per_multipole), so use 29.5 **/
/*****************************************************************************/
#define DIVISION_LIMBER_FULL 29.5

/**********************/ 
/** constant bias b0 **/
/**********************/
#define BIAS_CONST  2.0

/*******************************************************/ 
/** parameters to convert from E(B-V) to calibration  **/
/*******************************************************/ 
#define RV_DES_IBAND 1.595   // ratio of total to selective extinction, from Schlafly&Finkbeiner 2012
#define SZ_SLOPE 0.1 // s(z); faint-end slope of the mass function

/***********************/
/** z-binning of data **/
/***********************/
#define Z_BEGIN 0.0
#define DZ      0.20
#define Z_BINS    5 

#define Z_END  (Z_BEGIN+Z_BINS*DZ)

/***************************/
/** Survey specifications **/
/***************************/
#define FSKY         (5000.0 / 41253.0) 
#define OMEGA_SURVEY (FSKY * 4*M_PI) // solid angle

#define N_PER_SQARCMIN 20.0  // total n per square arcmin on sky (integrated over redshift)
#define N_PER_SR_TOT   (N_PER_SQARCMIN * 3600.0 * 41523.0 / (4*M_PI))  // total n per steradian on sky

/**************************/
/** to take derivatives ***/
/**************************/
#define FRAC_ERROR 0.02

#define ZMAX_D_TAB 3.0
#define DZ_D_TAB   0.0001
#define NSPLINE_REDSHIFT (lround(ZMAX_D_TAB/DZ_D_TAB)+1)

#define Q_NL 31.0 // used only at very low ell, so don't matter that we are using the inexact Tegmark et al... 
#define A_NL  1.4 // formula for nonlinear ps, since things are linear at those scales anyway!

/************************************************/
/*** Splining log(step) in k for P(k) and I(k) **/ 
/************************************************/
#define DLNK_PKSPLINE 0.1 // checked this is sufficient
#define DLNK_IKSPLINE 0.02 // checked this is sufficient

#define H0_hmpcinv  (1.0/2997.9) /* in h Mpc-1  */

/******************************************************************/
/** to converghe integrals, return jl'(x)=0 for x> (this value) **/
/******************************************************************/
#define MAX_ARG_RETURN_JLPRIME_ZERO 10.0

/***************/
/*** Global ****/ 
/***************/
double H0_mpcinv;
int  N_MASSFUN, N_GROWTH, Z_BIN_PAIRS;
/***************************/
/** jl'(x) and its spline **/
/***************************/
double  **xx_jlprime, **yy_jlprime, **y2_jlprime;

int WINDOW , NZ_CHOICE, NSPLINE_JLPRIME_NR;
double jl_prime(int l, double x);
void spline_jlprime_range_ell   (int lmax);
void spline_jlprime_range_ell_NR(int lmax);
double jlprime_tab   (int l, double x);
double jlprime_NR_tab(int l, double x);
gsl_spline *jlprime_spline[1000]; // to be safe!
double Wl_glob[101];         // maximum LMAX reasonably expected
double Wl_nhat_sn_glob[101][501]; // N_SN, N_LMAX
double MAX_ARG_JLPRIME_SPLINE, MIN_ARG_JLPRIME_SPLINE;
int COUNT_ABOVE, COUNT_IN, COUNT_BELOW;

/**************************/
/** for random variables **/
/**************************/

double ***alm_RE, ***alm_IM;

/**************************************/
/** For the Smith et al nonlinear PS **/
/*************************************/
double sigma_square(double z, double R);
void dfunc_sigsq(double x, double *f, double *fprime);
void initiate_smith(double z);
double KNL_SMITH, NEFF_SMITH, C_SMITH;
double ps_nonlin_smith(double k);

/************************************/
/** spline of (P(k) and T(k) stuff **/
/************************************/
int NSPLINE, NSPLINE_PK_CAMB, NSPLINE_JLPRIME_NR, NSPLINE_TK, NSPLINE_IK, N_L2, **l2_arr, *l2_1D_arr, *lm_index_arr, **l2m2_index_arr, CHOICE; 
gsl_spline *Pkspline, *Tkspline, *Pk_CAMB_spline, *r_spline, *W_spline[Z_BINS];
gsl_spline  *NL_LIN_CAMB_ratio_spline[Z_BINS], *NL_LIN_my_Smith_ratio_spline;
gsl_spline  *bias_spline[Z_BINS];
gsl_interp_accel *acc, *acc1, *acc2, *acc_Tk,  *acc_Pk_CAMB, *acc_bias, *acc3, *acc4, *acc5;
gsl_spline  *Ik_spline[Z_BINS], *Ik_spline_Tyson;
gsl_interp_accel *acc_Ik, *acc_Ik_Tyson;

int L_BINS, DIP_OR_QUAD, LMAX_CIJ_THETA, PRINT_FLAG;
double **Z_ARR, *n_per_sr, *dl_arr, *l_arr;

/*****************/
/**   for FoMs  **/
/*****************/
double FoM_DE(double **Finv);
double chisq_bias_w0wa(double *param_bias, double **Finv);

/********************************/
/* Kosowsky's hyperjl routines **/
/********************************/
extern double phi_recurs(int l, int K, double beta, double chi);
extern double phi_langer(int l, int K, double beta, double chi);

#ifdef __ultrix
extern double acosh(double);
#endif

double qintegral(double sin_K, double alpha, int K);
double airy_ai(double x);
double airy_aip(double x);
double polevl(double x, double *coef, int N);
double p1evl(double x, double *coef, int N);
void report_error(char error_text[]);

/*********************************/
/* For non-Limber power spectrum */
/*********************************/
struct I_params                  {double a; double b; int c;};
struct I_Tyson_params            {double a; double b;};
struct Cij_noLimber_params       {double a; double b; double c;};
struct Cij_noLimber_Tyson_params {double a;};
void spline_Ik     (double ell);
void spline_Ik_Tyson(double ell);
double I_tab      (int i, double k);
double I_Tyson_tab(double k); 
double I_integrand      (double z, void *p);
double I_Tyson_integrand(double z, void *p);
double I_int       (double ell, double k, int i);
double I_Tyson_int (double ell, double k);
double Cij_noLimber_integrand      (double log10k, void *p);
double Cij_noLimber_Tyson_integrand(double log10k, void *p);
double Cij_noLimber_int      (int do_nonlin, int i, int j);
double Cij_noLimber_Tyson_int(int do_nonlin, double Om, double w0, double wa);

/*********************/
/* for Cl */
/*********************/
double bias_simple(double fnl, double k, double b0);
void calculate_Cl(int firstpass,double om, double w0, double wa, double A, double n, double fnl,
                  double b0, double **Cl_ij_matrix, double **Cl_all_ell);
void tests_plots_Cl(int firstpass,
                    double omega_m, double w0, double wa, double A, double n, double fnl, double b0);

double Cl_vel_integrand(double log10k, void *p) ; 

double vel_integrand(double log10k, int ell, double z);

double Cij_theta_vel_integrand(double log10k, void *p);
double Cij_theta_gsl_int(int i, int j, double z1, double z2, double costh, double omega_m, double w0, double wa);


void run_camb_get_Tk_friendly_format(int do_nonlinear, double omega_m, double omhh_local, double obhh_local, double ns, 
                                     double dn_dlnk, double A_k_WMAP, double w0);
/* integration stuff */
struct Cij_params {int a; double b; double c;};
struct Cij_theta_params {int a; int b; double c; double d; double  e; double f; double g; double h;};

gsl_histogram *Tl_over_Cl_hist[L_MAX_CALIB];

/*********************************/
/** for SN and their covariance **/
/*********************************/
int  N_SN, NG_LIKE_CORRECTION, MARG_OVER_SCRIPTM;
void read_SN_pos(double **SN_z_th_phi, char *filename);
void read_pos_noise(int d, double **SN_z_th_phi, double *delta_m, double **Noise_Cov,
                  char *filename, char*file_noise_cov);
void order_SN_increasing_z(int nsn,double **SN_z_th_phi, double **Noise_C,double *delta_m);
void order_whole_SN_file_in_z(char *file_in, char *file_out);
void calculate_Cov_vel_of_SN(int nmax, double **SN_z_th_phi, double **Signal_SN, 
                             double omega_m, double w0, double wa);
void calculate_Signal_given_z_theta_arr(int ncosalpha, int nz, double *cosalpha_arr, double *z_arr, double ***Signal_tensor,
                                        double omega_m, double w0, double wa);
double getang_spher_cos(double theta1, double phi1, double theta2, double phi2);
void print_S_vs_z_fixed_theta(double z, double omega_m, double w0, double wa, char *filename) ; 
    
/***********************************/
/** for gals and their covariance **/
/***********************************/
void calculate_Cov_vel_of_gal(int nmax, double **SN_z_th_phi, double **Signal_SN,
                             double omega_m, double w0, double wa);
void read_gal_data(int nskip, int ngal, double **gal_z_th_phi, double *delta_m, double **Noise_Cov, char *filename);
void read_signal_matrix( double **Signal_SN, char *filename) ;
void order_whole_gal_file_in_z(int nskip, char *file_in, char *file_out);
/*********************************/
/** velocity  dispersion of SN  **/
/*********************************/
void   get_double_sky_integ_vs_ell(int ndir, int nside, int lmax, double **integrand_ell);
double vdisp_1D_integrand(double log10k, void *p);
double vdisp_1D_gsl_z1_z2_integrand(int ndir,double z1, double z2, double lmax, double omega_m, double w0, double wa);
double vdisp_1D(int ndir,  double zmin, double zmax, double lmax, double omega_m, double w0, double wa);

struct vdisp_1D_gsl_z1_z2_integrand_params {int a; double b; double c; double d; double e; double f; double g;};
void read_W_arrays(int lmax, int npix, char *root) ;

double *W_map_glob;
double **Wl_nhat_allsky_glob; // N_SN, N_LMAX
double **int_over_areas_vdisp_glob;

double sigmasq_vel_integrand(double log10k, void *p) ;
double sigmasq_vel(double R, double omega_m) ;
struct sigmasq_vel_integrand_params {double a; double b;};
int TOPHAT; // fiducial value is tophat, not exponential

void get_nz_of_JLA_SN(int sn, char *file);  // getting n(z) of actual SN
double nz_of_JLA(double z);  // call to n(z)
void print_histogram(gsl_histogram *h, char *filename);
gsl_histogram *nz_hist_glob;
/*********************************************/
/** geometric functions tabulation (for Cl) **/
/*********************************************/
double W_tab(int i, double z);
double Hsq(double z, double omega_m, double eqstate, double wprime);
double H_dimless(double z, double omega_m, double eqstate, double wprime);
double r_dimless_tab(double z), nz(double);
void spline_r(double omega_m, double eqstate, double wprime);
int file_eof_linecount(char *filename);
int count_columns_in_last_line_file(char *filename) ;
int nskip_first_line_of_file_or_not(char *filename) ;
void write_params_ini(FILE *ifp, int do_nonlinear, double H0, double obhh, double ochh, double ns, 
		      double dn_dlnk, double As, double w0);
/**************************************/
/** power spectrum tabulation and NL **/
/***************************************/
void spline_Pk(double z, double omega_m, double w0, double wa, double A, double n);
void spline_Tk_from_CAMB(double omega_m, double w0, double wa, double A, double n);
void spline_Pk_from_CAMB(double z, double omega_m, double w0, double wa, double A, double n);
void spline_Pk_from_CAMB_NR(double z, double omega_m, double w0, double wa, double A, double n);
void spline_bias(double fnl, double omega_m, double w0, double wa, double A, double n);

double NL_LIN_CAMB_ratio_tab(int i, double k), NL_LIN_my_Smith_ratio_tab(double k); 
void spline_my_Smith_NL_LIN_ratio_z(double z,double omega_m, double w0, double wa, double A, double n);
double Pk_tab(double k), Pk_CAMB_tab(double k), Pk_CAMB_NR_tab(double k);
double bias_tab(int i, double k);

double *xx_Pk , *yy_Pk, *y2_Pk;
/********************/
/* growth tabulated */
/********************/
double dlnD_dlna(double z_stop, double mat, double h, double w0, double w1);
double ff1(double a), ff2(double D, double Ddot, double lna, double w, double w1, double mat); 
double *xx, *yy , *y2, *xx_gr, *DD_gr, *aa_gr, *Dprime_gr, *gg_gr, *y2_gg, *y2_DD, *y2_Dprime;
double *dD_deta_gr, *y2_dD_deta;
double D_tab(double z), dD_da_tab(double z), g_tab(double z);
double dD_deta_tab(double z);
double D_comov(double zmax, double omega_m, double eqstate, double wprime);

/*********************/
/* for Fisher matrix */
/*********************/
void inverse(double **a, double **ai, int n);
double inverse_lndet_gsl_LU(double **A, double **Ainv, int n) ; 
void mat_mat_prod(double **a, double **b, double **c,int dim);
double trace(double **a, int dim);
double Fisher(double **Signal, double **Noise);
void print_matrix(double **a, int n);
void print_matrix_nonsquare(double **a, int m, int n);
void file_print_matrix(double **a, int n, char *file);
void file_read_matrix(double **a, int n, char *file) ; 
int ij_to_alpha(int i, int j);
void check_condition_number_matrix(int n, double **matrix);
double log_det_matrix(int n, double **matrix);

double chisq_cov(int n, double A, double *deltam, double **Signal, double **Noise) ; 
double chisq_cov_input_Cinv(int iter, int n, double *mu, double ***Cov_3D, double *logdetCov_arr, double ***Covinv_3D);
double chisq_in_mu_marg_over_scriptM(int n, double *v, double **Cinv, double Mmin, double Mmax, double dM);
double  calculate_norm_like_A(int nsn, double Amin, double Amax, double dA,
                              double sig_int_sq_max, double sig_int_sq_min, double dsig_int_sq,
                              double *deltam, double **Signal, double **Noise, char *filename);

void i_to_i1i2(int i, int *i1, int *i2, int i1max, int i2max) ; 

void find_hist_68_95_99_CL(gsl_histogram *h_orig, double *peak_min_max);
/*****************************************************************/
/** mass function stuff - not needed except for bias() function **/
/*****************************************************************/
double f_Warren(double s);
double dsigma_dM(double z, double M);
double sigma(double, double), sigma8(double);
double rho_m(double);
double *xx_RMS, *y2_RMS, *yy_RMS;


/*************************************/
/* global functions for Wayne's P(k)**/
/*************************************/
void set_global_variables_and_growth(double z, double omega_m, double w0, double wa, double A, double n, int calc_growth);

double TF_general (double kk), TF_BBKS(double kk);
double ps_linear(double kk);
double omm_z(double z, double mat, double w0, double w1);

/**************************/
/* more global variables **/
/**************************/
double omega_matter, omega_lambda, omega_hdm, w0_glob, wa_glob, fnl_glob,
  A_k0002,kk_fid,A_k0002_fid, omega_baryon,redshift, scalar_index;

/* The following are set in TFmdm_set_cosm() */
double growth_k0, growth_to_z0, hhubble,obhh,omega_curv,omega_lambda_z,omega_matter_z, omhh,onhh;


/* The following are set in TFmdm_onek_mpc() */
double	growth_to_z, growth_k0_new, growth_nu;

/* Finally, TFmdm_onek_mpc() and TFmdm_onek_hmpc() give their answers as */
double  growth_func;
