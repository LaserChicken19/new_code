#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"recipes/nrutil.h" 
#include"recipes/nr.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include "declarations.h" 

void spline_r(double omega_m, double eqstate, double wprime)
{
    // DIMENSIONLESS comoving distance r(z)
    int  i=0;
    double z;
    double x[NSPLINE_REDSHIFT], y[NSPLINE_REDSHIFT] ;

    acc1 = gsl_interp_accel_alloc ();
    r_spline = gsl_spline_alloc (gsl_interp_cspline, NSPLINE_REDSHIFT);

    double sum=0; 
    for(z=0; z<=ZMAX_D_TAB+SMALL; z+=DZ_D_TAB)
    {
        x[i] = z;
        y[i] = sum;
        sum+=DZ_D_TAB/H_dimless(z, omega_m, eqstate, wprime);
        i++;
    }

    gsl_spline_init (r_spline, x, y, NSPLINE_REDSHIFT);
}
void spline_Tk_from_CAMB(double omega_m, double w0, double wa, double A, double n)
/********************************************/
/** this is to tabulate the CAMB T(k, z=0) **/
/********************************************/
{
    int i;
    double k, Tk;
    FILE  *ifp;

    NSPLINE_TK = file_eof_linecount("my_transfer_out.dat");

    double x[NSPLINE_TK], y[NSPLINE_TK];

    acc_Tk   = gsl_interp_accel_alloc ();
    Tkspline = gsl_spline_alloc (gsl_interp_cspline, NSPLINE_TK);
   
    ifp=fopen("my_transfer_out.dat", "r");
    for (i=0; i<=NSPLINE_TK-1; i++)  // make sure you start at zero!
    {
        fscanf(ifp, "%lf %lf\n", &k, &Tk);
        x[i] = log10(k);
        y[i] = log10(Tk);
    }
    fclose(ifp);
    gsl_spline_init (Tkspline, x, y, NSPLINE_TK);

}
void spline_Pk_from_CAMB(double z, double omega_m, double w0, double wa, double A, double n)
/************************************************************************/
/** this is to tabulate the CAMB P(k, z=0) (note, NOT Delta^2(k)       **/
/************************************************************************/
{
    int i;
    double k, Pk_CAMB;
    FILE  *ifp;
    
    NSPLINE_PK_CAMB = file_eof_linecount("matterpower.dat");
    /****************************************************************************/
    /** Check if the first line starts with hashtag (ie. has variable names)   **/
    /** and skip that first line if that's the case                            **/
    /****************************************************************************/    
    int N_SKIP = nskip_first_line_of_file_or_not("matterpower.dat");
    printf("In matterpower (GSL tab), skipping %d lines\n", N_SKIP);    
    NSPLINE_PK_CAMB -= N_SKIP;

    /******************************/
    /** continue with tabulating **/
    /******************************/        
    double x[NSPLINE_PK_CAMB], y[NSPLINE_PK_CAMB];

    acc_Pk_CAMB    = gsl_interp_accel_alloc ();
    Pk_CAMB_spline = gsl_spline_alloc (gsl_interp_cspline, NSPLINE_PK_CAMB);
   
    ifp=fopen("matterpower.dat", "r");
    /*****************************/
    /** skip the unwanted lines **/
    /*****************************/
    int N_BUFFER=1e6;
    char unwanted_lines[N_BUFFER];    
    for(i=0; i < N_SKIP; i++)
        fgets(unwanted_lines, N_BUFFER, ifp);    
    
    for (i=0; i<=NSPLINE_PK_CAMB-1; i++)  // make sure you start at zero!
    {
        fscanf(ifp, "%lf %lf\n", &k, &Pk_CAMB);
        x[i] = log10(k);
        y[i] = log10(Pk_CAMB);
    }
    fclose(ifp);
    gsl_spline_init(Pk_CAMB_spline, x, y, NSPLINE_PK_CAMB);

}
void spline_Pk_from_CAMB_NR(double z, double omega_m, double w0, double wa, double A, double n)
/************************************************************************/
/** uses Numerical Recipes, not gsl, routine to tabulate               **/
/** this is to tabulate the CAMB P(k, z=0) (note, NOT Delta^2(k)       **/
/************************************************************************/
{
    int i;
    double k, Pk_CAMB;
    FILE  *ifp;
    
    NSPLINE_PK_CAMB = file_eof_linecount("matterpower.dat");
    /****************************************************************************/
    /** Check if the first line starts with hashtag (ie. has variable names)   **/
    /** and skip that first line if that's the case                            **/
    /****************************************************************************/    
    int N_SKIP = nskip_first_line_of_file_or_not("matterpower.dat");
    printf("In matterpower (NR tab), skipping %d lines\n", N_SKIP);    
    NSPLINE_PK_CAMB -= N_SKIP;

    /******************************/
    /** continue with tabulating **/
    /******************************/
    
    xx_Pk = dvector(1, NSPLINE_PK_CAMB);
    yy_Pk = dvector(1, NSPLINE_PK_CAMB);
    y2_Pk = dvector(1, NSPLINE_PK_CAMB);

    ifp=fopen("matterpower.dat", "r");
    /*****************************/
    /** skip the unwanted lines **/
    /*****************************/
    int N_BUFFER=1e6;
    char unwanted_lines[N_BUFFER];    
    for(i=0; i < N_SKIP; i++)
        fgets(unwanted_lines, N_BUFFER, ifp);    

    for (i=1; i<=NSPLINE_PK_CAMB; i++)  // make sure you start at zero!
    {
        fscanf(ifp, "%lf %lf\n", &k, &Pk_CAMB);
        xx_Pk[i] = log10(k);
        yy_Pk[i] = log10(Pk_CAMB);
    }
    fclose(ifp);

    my_spline(xx_Pk, yy_Pk, NSPLINE_PK_CAMB, 0, 1e31, y2_Pk);
}
    
double Pk_CAMB_NR_tab(double k)
{
    double logPk;
  
    my_splint(xx_Pk, yy_Pk, y2_Pk, NSPLINE_PK_CAMB, log10(k), &logPk);
   return pow(10, logPk);
}
void spline_Pk(double z, double omega_m, double w0, double wa, double A, double n)
{
    int  i=0;
    double x[NSPLINE], y[NSPLINE], k;
    
    acc = gsl_interp_accel_alloc ();
    Pkspline = gsl_spline_alloc (gsl_interp_cspline, NSPLINE);

    /****************************************************************************/
    /** first, initiate the global cosmological parameters and growth function **/
    /****************************************************************************/
    set_global_variables_and_growth(z, omega_m, w0, wa, A, n, 1);
    
    /**************************/
    /** Spline P(k) as usual **/
    /**************************/    
    for(k=K_START_INTERP; k<=K_END_INTERP; k+=DLNK_PKSPLINE*k)
    {
        x[i] = log10(k);
        y[i] = log10(ps_linear(k)); 
        i++;
    }
    gsl_spline_init (Pkspline, x, y, NSPLINE);
}


void dfunc_sigsq(double x, double *f, double *fprime)
{
    double dummy, dummy1, dx;
  
    dx=0.02*x;    
    dummy =  sigma_square(redshift, x) - 1.0;
    dummy1= (sigma_square(redshift, x+dx) - sigma_square(redshift, x))/dx;
    
    *f=dummy;
    *fprime=dummy1;
}

void initiate_smith(double z)
{
    double dR, r, R, d1, d2, sigsq;

  redshift = z; // goes into sigma_square just below
  
  R=my_rtsafe(dfunc_sigsq, 0.1, 20.0, 0.01); /* in Mpc; find R for which sigma^2(z, R)=1.0 */
  dR=0.02*R;

  r = R * hhubble;              /* in Mpc/h */
  KNL_SMITH=1.0/r;              /* in h/Mpc */
  
  sigsq=sigma_square(z, R); // recall, need R in sigma_square
  d1=(sigma_square(z, R+dR)-sigma_square(z, R-dR))/2.0/dR;
  d2=(sigma_square(z, R+dR)+sigma_square(z, R-dR)-2*sigma_square(z, R))/(dR*dR);

  NEFF_SMITH = -R/sigsq*d1-3.0;                 // double checked in Dec 2011

  C_SMITH = -R*(d1/sigsq + R/sigsq*d2 - R/(sigsq*sigsq)*d1*d1); // double checked in Dec 2011

}

void spline_my_Smith_NL_LIN_ratio_z(double z, double omega_m, double w0, double wa, double A, double n)
{
    /***********************************************************************************/
    /** splines Smith et al the NL/LIN ratio at all redshift and all k, using MY code **/
    /***********************************************************************************/

    int  i ;
    double k;

    double x[NSPLINE], y[NSPLINE];    
    acc3 = gsl_interp_accel_alloc ();

    NL_LIN_my_Smith_ratio_spline = gsl_spline_alloc (gsl_interp_cspline, NSPLINE);

    /*************************************************************************************************/
    /** initiate the global cosmological parameters - although they are already here from Pk_spline **/
    /** save time by NOT computing thed growth function  - already there from Pk_spline             **/
    /** but you do need growth function EVALUATED at this z, so run this                            **/
    /*************************************************************************************************/
    set_global_variables_and_growth(z, omega_m, w0, wa, A, n, 0);

    initiate_smith(z);
        
    i=0;
    for(k=K_START_INTERP; k<=K_END_INTERP; k+=DLNK_PKSPLINE*k)
    {
        x[i] = log10(k);
        y[i] = ps_nonlin_smith(k)/ps_linear(k);
        i++;
    }

    gsl_spline_init (NL_LIN_my_Smith_ratio_spline, x, y, NSPLINE);

}
void spline_bias(double fnl, double omega_m, double w0, double wa, double A, double n)
{
    /***********************************************************************************/
    /** splines Smith et al the NL/LIN ratio at all redshift and all k, using MY code **/
    /***********************************************************************************/

    int  i, bin;
    double z, k;

    double x[NSPLINE], y[NSPLINE];    
    acc_bias = gsl_interp_accel_alloc ();

    for(bin=1; bin<=Z_BINS; bin++)
    {        
        bias_spline[bin-1] = gsl_spline_alloc (gsl_interp_cspline, NSPLINE);

        z = 0.5*(Z_ARR[bin][1]+Z_ARR[bin][2]);

        /*************************************************************************************************/
        /** initiate the global cosmological parameters - although they are already here from Pk_spline **/
        /** save time by NOT computing thed growth function  - already there from Pk_spline             **/
        /** but you do need growth function EVALUATED at this z, so run this                            **/
        /*************************************************************************************************/
        set_global_variables_and_growth(z, omega_m, w0, wa, A, n, 0);

        i=0;
        for(k=K_START_INTERP; k<=K_END_INTERP; k+=DLNK_PKSPLINE*k)
        {
            x[i] = log10(k);
            y[i] = bias_simple(fnl, k, BIAS_CONST);
            i++;
        }

        gsl_spline_init (bias_spline[bin-1], x, y, NSPLINE);
    }

}
double ps_nonlin_smith(double k)
/* k is in h/Mpc, and so should be KNL_SMITH */
{
  double a, b, c, alpha, beta, gamma, mu, nu, f_1b, f_2b, f_3b;
  double n, C, lin, Delta_H, Delta_Q, y;

  n=NEFF_SMITH;
  C=C_SMITH;

  a=pow(10.0, 1.4861 + 1.8369*n + 1.6762*n*n + 0.7940*n*n*n +
	0.1670*n*n*n*n - 0.6206*C);

  b=pow(10.0, 0.9463 + 0.9466*n + 0.3084*n*n - 0.9400*C);
  c=pow(10.0, -0.2807 + 0.6669*n + 0.3214*n*n - 0.0793*C);

  gamma=0.8649 + 0.2989*n + 0.1631*C;
  alpha=1.3884 + 0.3700*n - 0.1452*n*n;
  beta =0.8291 + 0.9854*n + 0.3401*n*n;

  mu=pow(10.0, -3.5442 + 0.1908*n);
  nu=pow(10.0, 0.9589 + 1.2857*n);

  f_1b=pow(omega_matter_z, -0.0307);
  f_2b=pow(omega_matter_z, -0.0585);
  f_3b=pow(omega_matter_z,  0.0743);

  lin=ps_linear(k);
  y=k/KNL_SMITH;

  Delta_Q=lin*pow(1+lin, beta)/(1+alpha*lin)*exp(-(y/4+y*y/8));
  
  Delta_H=a*pow(y, 3*f_1b)/(1 + b*pow(y, f_2b) + pow(c*f_3b*y, 3-gamma));

  Delta_H = Delta_H / (1 + mu/y + nu/(y*y));

  return(Delta_Q + Delta_H);
}
double Pk_tab(double k) /* really, Delta^2(k) */
{
    double logPk = gsl_spline_eval (Pkspline, log10(k), acc);

    return pow(10, logPk);
}
double Pk_CAMB_tab(double k) /* Note, this is the traditional P(k), NOT Delta^2 */
{
    if (k > K_END_INTERP)
        printf("OVER bound of interpolated k!  k=%f\n", k);
    double logPk_CAMB = gsl_spline_eval (Pk_CAMB_spline, log10(k), acc_Pk_CAMB);

    return pow(10, logPk_CAMB);
}
double NL_LIN_CAMB_ratio_tab(int i, double k) 
{
    double ratio = gsl_spline_eval (NL_LIN_CAMB_ratio_spline[i-1], log10(k), acc4);

    return ratio;
}
double NL_LIN_my_Smith_ratio_tab(double k) 
{
    double ratio = gsl_spline_eval (NL_LIN_my_Smith_ratio_spline, log10(k), acc3);

    return ratio;
}
double bias_tab(int i, double k) 
{
    double b = gsl_spline_eval (bias_spline[i-1], log10(k), acc_bias);

    return b;
}
double r_dimless_tab(double z) //r(z)
{
    double rz = gsl_spline_eval (r_spline, z, acc1);

    return rz;
}
double W_tab(int i, double z) //w_i(z)
{
    return gsl_spline_eval (W_spline[i-1], z, acc2);
}

double Cl_vel_integrand(double log10k, void *p)
{ 
    // velocity power spectrum, following the note of Hui
    double k;

    struct Cij_params * params = (struct Cij_params *)p;

    int ell   = (params->a);
    double z1 = (params->b);
    double z2 = (params->c);

    k = pow(10.0, log10k);  // units are h Mpc^{-1}

    double chi1 = r_dimless_tab(z1)/H0_hmpcinv;
    double chi2 = r_dimless_tab(z2)/H0_hmpcinv;
    
    return( 2/M_PI *
            log(10.0) * k *                         // conversion: they have dk, I want dlog10(k)
            2*M_PI*M_PI * Pk_tab(k) * pow(k, -3) *  // P(k)
            //Pk_CAMB_tab(k) *                        // NLPk from CAMB, leads to small diff @ higher ell, Apr 2015
            dD_deta_tab(z1) * dD_deta_tab(z2)     * // units are [h Mpc^{-1}]^2
            jl_prime(ell, k*chi1) *  jl_prime(ell, k*chi2)
            );
}

double jl_prime(int l, double x)
{
    if (l == 0) return (-gsl_sf_bessel_jl(1, x));
    else 
        return (gsl_sf_bessel_jl(l-1, x) - (l+1)*gsl_sf_bessel_jl(l, x)/x);
}
void spline_jlprime_range_ell(int lmax)
{
    // jl'(x)
    int  i=0, l, NUM_STEPS=10000;
    double arg;
    double x[NUM_STEPS+1], y[NUM_STEPS+1];
    
    double  min_arg=1.0e-3, max_arg=1000.0;

    double dln_arg = 1.0/NUM_STEPS*log(max_arg/min_arg);

    MIN_ARG_JLPRIME_SPLINE=min_arg;
    MAX_ARG_JLPRIME_SPLINE=max_arg;

    acc2 = gsl_interp_accel_alloc ();

    for(l=0; l<=lmax; l++)
    {
        jlprime_spline[l] = gsl_spline_alloc (gsl_interp_cspline, NUM_STEPS+1);

        i=0;

        for(i=0; i<=NUM_STEPS; i++)
        {
            arg = min_arg * exp(i*dln_arg);        
            x[i] = arg;
            if (arg < 0.5*l)   // at arg<l, jl'(arg) hugs zero like cra-a-zy and GSL can give an underflow (Oct 2016)
                y[i] = 0;
            else
                y[i] = jl_prime(l, arg);
        }
        gsl_spline_init (jlprime_spline[l], x, y, NUM_STEPS+1);
    }// ends for that l

}
double jlprime_tab(int l, double x)
{
    if (x < MIN_ARG_JLPRIME_SPLINE) COUNT_BELOW++;
    else if (x > MAX_ARG_JLPRIME_SPLINE) COUNT_ABOVE++;
    else COUNT_IN++;

    //printf("l=%d, x=%f; min=%f max=%f\n", l, x, MIN_ARG_JLPRIME_SPLINE, MAX_ARG_JLPRIME_SPLINE);
    
    if ((x < MIN_ARG_JLPRIME_SPLINE) || (x > MAX_ARG_JLPRIME_SPLINE))
        return jl_prime(l, x);
    else        
        return gsl_spline_eval (jlprime_spline[l], x, acc2);
}
void spline_jlprime_range_ell_NR(int lmax)
{
  /************************************************************/
  /**                     jl'(x)                             **/
  /** uses Numerical Recipes, not gsl, routine to tabulate   **/
  /** note that you want xx[num_ell][num_interp] to feed in  **/
  /**           see code test_NR_2d_spline                   **/
  /************************************************************/

    int  i, l;

    NSPLINE_JLPRIME_NR=10000; // defined in declarations
    double arg;

    /******************************************************/
    /** arrays of tabulation arrays for each ell         **/
    /** key part: first count in ell, then in argument   **/
    /******************************************************/
    xx_jlprime = dmatrix(1, lmax+1,1, NSPLINE_JLPRIME_NR);
    yy_jlprime = dmatrix(1, lmax+1,1, NSPLINE_JLPRIME_NR);
    y2_jlprime = dmatrix(1, lmax+1,1, NSPLINE_JLPRIME_NR);
    
    double  min_arg=1.0e-3, max_arg=1000.0;

    double dln_arg = 1.0/(NSPLINE_JLPRIME_NR-1)*log(max_arg/min_arg);

    MIN_ARG_JLPRIME_SPLINE=min_arg;
    MAX_ARG_JLPRIME_SPLINE=max_arg;

    for(l=0; l<=lmax; l++)
    {

        for(i=1; i<=NSPLINE_JLPRIME_NR; i++)
        {
            arg = min_arg * exp((i-1)*dln_arg);

            //printf("HERE %e (l=%d)\n", arg, l);
            xx_jlprime[l+1][i] = arg;
            if (arg < 0.1*l)   // at arg<<l, jl'(arg) hugs zero like cra-a-zy and GSL can give an underflow (Oct 2016)
                yy_jlprime[l+1][i] = 0;
            else
                yy_jlprime[l+1][i] = jl_prime(l, arg);
        }
        
        // spline each value of l
        my_spline(xx_jlprime[l+1], yy_jlprime[l+1], NSPLINE_JLPRIME_NR, 0, 1e31, y2_jlprime[l+1]);

    }// ends for that l

}
double jlprime_NR_tab(int l, double x)
{
    double jlp;

    if ((x < MIN_ARG_JLPRIME_SPLINE) || (x > MAX_ARG_JLPRIME_SPLINE))
        jlp = jl_prime(l, x);

    else
        my_splint(xx_jlprime[l+1], yy_jlprime[l+1], y2_jlprime[l+1], NSPLINE_JLPRIME_NR, x, &jlp);

return jlp;
}    

      
double vel_integrand(double log10k, int ell, double z)
{
    double z1=z;
    double z2=z;
    double k = pow(10.0, log10k); // proportional to H0

    double chi1 = r_dimless_tab(z1)/H0_hmpcinv;
    double chi2 = r_dimless_tab(z2)/H0_hmpcinv;
    
    return( 2/M_PI *
            log(10.0) * k *                        // conversion: they have dk, I want dlog10(k)
            2*M_PI*M_PI * Pk_tab(k) * pow(k, -3) * // P(k)
            dD_deta_tab(z1) * dD_deta_tab(z2)     *  // proportional to H0^2
            jl_prime(ell, k*chi1) *  jl_prime(ell, k*chi2)
            );
}    
double Cij_theta_vel_integrand(double log10k, void *p)
{ 
    // velocity power spectrum, following the note of Hui
    int ell;
    double k;
    double sum=0, fac;

    struct Cij_theta_params * params = (struct Cij_theta_params *)p;

    int     i = (params->a);
    int     j = (params->b);
    double z1 = (params->c);
    double z2 = (params->d);
    double costh = (params->e);
    double om = (params->f);
    double w0 = (params->g);       
    double wa = (params->h);

    k    = pow(10.0, log10k);  // units are h Mpc^{-1}
    
    double chi1 = r_dimless_tab(z1)/H0_hmpcinv;
    double chi2 = r_dimless_tab(z2)/H0_hmpcinv;

    double prefac0 = pow(5/log(10.0), 2);
    // checked very small diff in results if the starting 1s are zeroed out
    double prefac1 = 1 - (1+z1)/H_dimless(z1, om, w0, wa)/r_dimless_tab(z1);
    double prefac2 = 1 - (1+z2)/H_dimless(z2, om, w0, wa)/r_dimless_tab(z2);
    if (PRINT_FLAG == 1) printf("%f %f %d\n",  prefac0, prefac1, NSPLINE_PK_CAMB);
    double prefac3 =  1/(2*M_PI*M_PI) *
        log(10.0) * k *                        // conversion: they have dk, I want dlog10(k)
        Pk_CAMB_NR_tab(k) *                       // P(k) from CAMB
        dD_deta_tab(z1) * dD_deta_tab(z2)      // units are [h Mpc^{-1}]^2
        ;

    int lmax;

    /////////////////////////////////////////////////////////////////////////
    // 16 October 2016: found sufficient conditions on lmax;
    // drastically smaller lmax in all cases except aligned, high-z objects
    /////////////////////////////////////////////////////////////////////////
    //if ((z1 > 0.05) && (z2 > 0.05) && (costh > 0.95)) lmax = LMAX_CIJ_THETA;
    if (costh > 0.95)  lmax = LMAX_CIJ_THETA;
    else               lmax = 20.0;
    //FIX
    //lmax=10;
    
    for (ell=0; ell<=lmax; ell++)
    {
        if (WINDOW == 3) // no mean of magnitude subtracted
            fac = gsl_sf_legendre_Pl(ell,costh);
        if (WINDOW == 2) // mean subtracted with arbitrary window
        {
            fac = gsl_sf_legendre_Pl(ell,costh) -
                4*M_PI/(2*ell+1) * Wl_nhat_sn_glob[ell][i] -
                4*M_PI/(2*ell+1) * Wl_nhat_sn_glob[ell][j] +
                4*M_PI           * Wl_glob[ell];
        }
        if (WINDOW == 1) // mean subtracted with window = all sky
        {
            if (ell == 0) // W_l and W_l(nhat) are nonzero only for ell=0
                fac = gsl_sf_legendre_Pl(ell,costh) - 2 + 1;
            else 
                fac = gsl_sf_legendre_Pl(ell,costh);
        }

        //sum += (2*ell+1) * jl_prime(ell, k*chi1) *  jl_prime(ell, k*chi2) * fac;
        // checked same result with jlprime_NR_tab as with jlprime_tab and jl_prime - Apr 2016
        sum += (2*ell+1) * jlprime_NR_tab(ell, k*chi1) *  jlprime_NR_tab(ell, k*chi2) * fac;
    }

    return(prefac0*prefac1*prefac2*prefac3*sum );
}
double Cij_theta_gsl_int(int i, int j, double z1, double z2, double costh, double omega_m, double w0, double wa)
{
    /*********************************/
    /** Following the Hui_note.pdf  **/
    /*********************************/

    /*******************************************************************/
    /**             !!!!  I M P O R T A N T !!!!                      **/
    /** Turn off the error handler since things are already accurate  **/
    /*******************************************************************/
    gsl_set_error_handler_off ();
    
    /* for gsl integration */
    gsl_integration_workspace * work = gsl_integration_workspace_alloc (1000);
    double result, error;

    gsl_function F;
    F.function = &Cij_theta_vel_integrand;
           
    struct Cij_theta_params params = {i, j, z1, z2, costh, omega_m, w0, wa};
    F.params = &params;

    //printf("%f %f\n", LOG10_KMIN_CIJ_INT, LOG10_KMAX_CIJ_INT);
    /* run the integrator! */
    /************************************************************************************/
    /*** Oct 16, 2016: found that need high kmax only for very low z, and vice versa.  **/
    /************************************************************************************/
    double log10_kmax_func_of_z;
    if     ((z1 < 0.005) && (z2 < 0.005))  log10_kmax_func_of_z = 2.0;
    else if ((z1 < 0.03) && (z2 < 0.03))   log10_kmax_func_of_z = 1.0;
    else                                   log10_kmax_func_of_z = 0.5; // agrees with LOG10_KMAX_CIJ_INT
    //FIX
    //log10_kmax_func_of_z=0.0;
    
    gsl_integration_qag (&F, LOG10_KMIN_CIJ_INT, log10_kmax_func_of_z,
                         1.0e-7, 1.0e-2, 1000, 1, work, &result, &error);  // -2 is good enough for rel, given LMAX and KMAX unc
    //printf ("result          = % .18f\n", result);
    //printf ("estimated error = % .18f\n", error);
    //printf ("intervals =  %d\n", work->size);
  
    gsl_integration_workspace_free (work);

    return(result);
}

