#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"recipes/nrutil.h" 
#include"recipes/nr.h"
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include "declarations.h"

double D_tab(double z)
{
  double gro;
  
   my_splint(xx_gr, DD_gr, y2_DD, N_GROWTH,z,&gro);
   return(gro) ;
}
double dD_da_tab(double z)
{
    //////////////
    // !! dD/da !!
    //////////////
    double dD_da;
  
    my_splint(xx_gr, Dprime_gr, y2_Dprime, N_GROWTH, z, &dD_da);
    return(dD_da) ;
}
double dD_deta_tab(double z)
{
    ///////////////////////////
    /// has the same unit as H
    ///////////////////////////
    double dD_deta;
  
    my_splint(xx_gr, dD_deta_gr, y2_dD_deta, N_GROWTH, z, &dD_deta);
    return(dD_deta) ;
    
}
    
double g_tab(double z)
{
  double gro;
  
   my_splint(xx_gr, gg_gr, y2_gg, N_GROWTH,z,&gro);
   return(gro) ;
}
double dsigma_dM(double z, double M)
{
    /************************/
    /*** 2-sided dsigma/dM  */
    /************************/

    double sig1, sig2, dM=0.1*M;

    my_splint(xx_RMS, yy_RMS, y2_RMS, N_MASSFUN, M-dM, &sig1);
    sig1=sig1*D_tab(z);
    
    my_splint(xx_RMS, yy_RMS, y2_RMS, N_MASSFUN, M+dM, &sig2);
    sig2=sig2*D_tab(z);

    return((sig2-sig1)/(2*dM));
}
double sigma_square(double z, double R)
/**************************************************************/
/** with GAUSSIAN filtering, as per AppendixC of Smith et al **/
/** R is accepted in Mpc,and converted to r in Mpc/hh  below **/
/** called in initiate_smith, before which globals have been set already **/
/*************************************************************************/
{  
    double r,  sigma_sq, k, dk;

    r = R * hhubble;  /* covert to Mpc/h:- need it since k below is in h/Mpc */
  
    sigma_sq=0;
  for (k=0.001*(1.0/R); k<=100.0*(1.0/R); k+=dk)
    {
      /*printf("%e %e\n", k, R);*/
      sigma_sq+=dk/k*ps_linear(k)*exp(-k*k*r*r);
      dk=0.01*k;   
    }

  // fix on Sep 27, 2016: you had been forgetting to mu;tiply by D(z)!!!
  return(D_tab(z)*D_tab(z)*sigma_sq);
}
double sigma8(double z)
/*
   k is in h Mpc^{-1}  */
{
  double r, sigma_sq, k, dk;
  float spher_bessel, dummy, dummy1, dummy2;

  r=8.0;  /* in h^{-1} Mpc */

  dk=0.0001;  
  sigma_sq=0;
  for (k=0.001; k<=100.; k+=dk) // 100 is sufficient for infinite accuracy - 31 Oct 2011
    {
      sphbes(1, k*r, &spher_bessel, &dummy, &dummy1, &dummy2);
      sigma_sq+=dk/k*ps_linear(k)*pow(3*spher_bessel/(k*r), 2);
      /*
	printf("%e %f %e %f %f %f %f %f %e\n", k, growth_to_z0, 
	       ps_linear(k), w, hhubble, k*hhubble/H0_mpcinv,
	       omega_matter_z, omega_lambda_z, TF_general(k));
      */
      dk=0.01*k;   /* 0.01 -> result accurate to 0.3% */
    }
  /*
  printf( "Sigma is %f\n", sqrt(sigma_sq));
  */

  // fix on Sep 27, 2016: you had been forgetting to mu;tiply by D(z)!!!
  return(D_tab(z)*sqrt(sigma_sq));
}
double rho_m(double z)
{
  return(CRIT_DENS*omega_matter*pow(1+z, 3));
}
double bias_simple(double fnl, double k, double b0)
{
    /*****************************************************************/
    /*              Used in splining b(fnl, k)                       */
    /*             takes Gaussian  b0 as input parameter             */
    /* k is given in h Mpc^{-1} and only affects the fnl \neq 0 case */
    /* Requires pre-computation of transfer function and setting of global cosmological parameters **/
    /*****************************************************************/
    double Deltab;

    Deltab = (b0-1) * fnl *
        (3*omega_matter*H0_hmpcinv*H0_hmpcinv*DELTA_C) /
        (growth_k0_new * TF_general(hhubble*k) * k * k);

    return b0 + Deltab;
}
                                            
double ps_linear(double kk)
/* Returns LINEAR Delta^2(k)
   with transfer function of Eisenstein & Hu (1997)  */
/* this is same as 1/(2PI^2)*k^3P(k) */
/* kk is given in h Mpc^{-1}*/
/* NOTE: every time this function is called, global variables (and z) have already been set */
/* Returns the density-weighted CDM+Baryon power spectrum */
/* transfer function is given in TF_general */
{ 
  double  transfer, k_fid;

  /* using Hu & Einsenstein */
  kk=kk*hhubble; /* revert to "ordinary" kk in Mpc^{-1} */

  transfer=TF_general(kk);

  k_fid=0.002;    /* WMAP expands around this; in Mpc^{-1} */

  /********************************************************************/
  /* using the normalization norm=A_k0002 defined at k=0.002 Mpc^{-1} **/
  /* note power spectrum \propto A, not A^2                           **/
  /*********************************************************************/

  return(4.0/25*A_k0002*pow(omega_matter, -2)*pow(kk/k_fid, scalar_index-1)*
	 pow(kk/H0_mpcinv, 4)*pow(growth_k0_new*transfer, 2));

} 

double TF_general(double kk)
{
  /* kk is given in Mpc^{-1} */
  double transfer;

  // Using Hu-Eisenstein
  //transfer=TFmdm_onek_mpc(kk);

  /* takes kk in h*Mpc^{-1} */
  transfer = pow(10.0, gsl_spline_eval (Tkspline, log10(kk/hhubble), acc_Tk));
  
  /**************************/ 
  /****   using BBKS   ******/
  //transfer =TF_BBKS(kk); 
  /**************************/ 
  
  return transfer;
}
void set_global_variables_and_growth(double z, double omega_m, double w0, double wa, double A, double n, int calc_growth)
{

    redshift = z;
    A_k0002 = A;
    scalar_index = n;
    
    omega_matter = omega_m;    
    omega_lambda = 1.0-omega_matter;    

    hhubble      = sqrt(omhh/omega_m);
    omega_baryon = obhh*pow(hhubble, -2);
    omega_curv   =  0.0;

    /************************************************/
    /** first, integrate to get the linear growth! **/
    /************************************************/
    if (calc_growth == 1)
        dlnD_dlna(100.0, omega_m, hhubble, w0, wa); 
    
    /************************************************************/
    /** From the direct integraton of the growth 2nd order eq: **/
    /************************************************************/
    growth_to_z0  = D_tab(redshift);
    growth_k0     = g_tab(redshift); /* this is absolute growth starting at z=Z_EQ */

    /********************************************************/
    /** this is really a*g(a); a at high z and 0.75 at z=0 **/ 
    /********************************************************/
    growth_k0_new = growth_k0/Z_EQ;

    
    /********************************************************/
    /** this is really g(a); one at high z and 0.75 at z=0 **/ 
    /********************************************************/
    growth_func   = growth_k0_new*(1.0+redshift); 

    /*******************************/
    /** Omega_M(z) for (w0, wa)   **/ 
    /*******************************/
    omega_matter_z=
        omega_matter*pow(1+redshift, 3)/
        (omega_matter*pow(1+redshift, 3)+omega_curv*pow(1+redshift, 2)+
         omega_lambda*pow(1+redshift, 3*(1+w0+wa))*
         exp(-3*redshift*wa/(1+redshift)));

}

double TF_BBKS(double kk)
{  /* kk is given in Mpc^{-1}!*/
  /* the transfer function using bbks formulae */

  double qval, Gamma, om, ob;

  // 11 Feb 2016: ob and omega_matter defined in galcov_main
  
  om=omega_matter;
  ob=omega_baryon;
  //hhubble=sqrt(omhh/om);

  //ob=obhh/hhubble/hhubble;
  
  Gamma=om*hhubble*exp(-ob-1.3*ob/om);
  qval = kk/(Gamma*hhubble);

  return (log(1+2.34*qval)/(2.34*qval)*
    pow(1 + 3.89*qval + pow(16.1*qval,2.0) + pow(5.46*qval,3.0) +
        pow(6.71*qval,4.0), -1.0/4) );
}
double ff1(double a)
{
  return(a);
}

double ff2(double D, double Ddot, double lna, double w, double w1,
          double mat) 
{
  double a, z;
  a=exp(lna);
  z=1.0/a-1;

  return(-Ddot*Ddot - Ddot*(0.5-1.5*(w+w1*z/(1+z))*(1-omm_z(z, mat, w, w1))) + 
	 1.5*omm_z(z, mat, w, w1));
}

double omm_z(double z, double mat, double w0, double w1)
{

  return(mat*pow(1+z, 3)/(mat*pow(1+z, 3) +
			  (1-mat)*pow(1+z, 3*(1+w0+w1))*exp(-3*w1*z/(1+z))));
}
double dlnD_dlna(double z_stop, double mat, double h, double w0, double w1)
{
  /***************************************************/
  /** integrates the second order equation in ln(D) **/
  /** returns the value of dlnD/dlna *****************/
  /***************************************************/
  double A_STOP,  A_INIT;
  double lna, del, deldot, deldot_return, dlna,k1_1,k1_2,k1_3,k1_4,k2_1, k2_2;
  double ddel, ddeldot, k2_3, k2_4, a, D, z, Hdim;
  int i;


  del=0.0; /* so that starting D=exp(del)=1; important only to get g(z)  */
  /* IMPORTANT-with reasonable guess, will lock into desired sol. faster */  
  deldot=1.0; 

  A_INIT=1.0/Z_EQ;   // doesn't matter what exactly Z_EQ is, as long as it's high z
  A_STOP=1.0/(1+z_stop);

  dlna=0.001; 

  i=0;
  for(lna=log(A_STOP); lna<=0; lna+=dlna) i++;
  N_GROWTH=i;

  //printf("N_GROWTH=%d\n", N_GROWTH);

  xx_gr=dvector(1, N_GROWTH);
  DD_gr=dvector(1, N_GROWTH);
  gg_gr=dvector(1, N_GROWTH);
  y2_gg=dvector(1, N_GROWTH);
  y2_DD=dvector(1, N_GROWTH);

  aa_gr=dvector(1, N_GROWTH);
  Dprime_gr=dvector(1, N_GROWTH);
  y2_Dprime=dvector(1, N_GROWTH);

  dD_deta_gr=dvector(1, N_GROWTH);
  y2_dD_deta=dvector(1, N_GROWTH);

  for(lna=log(A_INIT); lna<=log(A_STOP); lna+=dlna) 
    {
      
      k1_1=dlna*ff1(deldot) ;
      k2_1=dlna*ff2(del, deldot, lna, w0, w1, mat);
      
      k1_2=dlna*ff1(deldot + k2_1/2);
      k2_2=dlna*ff2(del + k1_1/2, deldot+k2_1/2, lna+dlna/2, w0, w1, mat);
      
      k1_2=dlna*ff1(deldot + k2_1/2);
      k2_2=dlna*ff2(del + k1_1/2, deldot+k2_1/2, lna+dlna/2, w0, w1, mat);
      
      k1_3=dlna*ff1(deldot + k2_2/2);      
      k2_3=dlna*ff2(del + k1_2/2, deldot+k2_2/2, lna+dlna/2, w0, w1, mat);
      
      k1_4=dlna*ff1(deldot +k2_3);
      k2_4=dlna*ff2(del + k1_3, deldot + k2_3, lna+dlna,   w0, w1, mat);
      
      ddel    = (k1_1+2*k1_2+2*k1_3+k1_4)/6;
      ddeldot = (k2_1+2*k2_2+2*k2_3+k2_4)/6;
      del     = del + ddel;
      deldot  = deldot + ddeldot;
      
      a=exp(lna);
      D=exp(del);
    }
  deldot_return=deldot;
  
  i=N_GROWTH;
  /* contunue to present time in order to normalize D */
  for(lna=log(A_STOP); lna<=0; lna+=dlna) 
        {

          k1_1=dlna*ff1(deldot) ;
          k2_1=dlna*ff2(del, deldot, lna, w0, w1, mat);
          
          k1_2=dlna*ff1(deldot + k2_1/2);
          k2_2=dlna*ff2(del + k1_1/2, deldot+k2_1/2, lna+dlna/2, w0, w1, mat);
          
          k1_2=dlna*ff1(deldot + k2_1/2);
          k2_2=dlna*ff2(del + k1_1/2, deldot+k2_1/2, lna+dlna/2, w0, w1, mat);
          
          k1_3=dlna*ff1(deldot + k2_2/2);      
          k2_3=dlna*ff2(del + k1_2/2, deldot+k2_2/2, lna+dlna/2, w0, w1, mat);
          
          k1_4=dlna*ff1(deldot +k2_3);
          k2_4=dlna*ff2(del + k1_3, deldot + k2_3, lna+dlna,     w0, w1, mat);
          
          ddel    = (k1_1+2*k1_2+2*k1_3+k1_4)/6;
          ddeldot = (k2_1+2*k2_2+2*k2_3+k2_4)/6;
          del     = del + ddel;
          deldot  = deldot + ddeldot;

	  a=exp(lna);
	  D=exp(del);
	  
	  xx_gr[i]=1.0/a-1;
          aa_gr[i]=a;
	  gg_gr[i]=D;
	  /**printf("in dlnD %f %e\n", xx_gr[i], gg_gr[i]);*/
	  i--;
        }

  for(i=N_GROWTH; i>=1; i--)
      DD_gr[i]=gg_gr[i]/D;       /* normalized D array */
  for(i=N_GROWTH; i>=2; i--)
  {
      Dprime_gr[i-1]=(DD_gr[i-1]-DD_gr[i]) / (aa_gr[i-1]-aa_gr[i]);       /* dD/da */

      // get dD/deta; April 2015
      a = 0.5*(aa_gr[i-1]+aa_gr[i]);
      z = 1.0/a-1;
      Hdim = H_dimless(z,  mat, w0, w1)*H0_hmpcinv;                /* units: h Mpc^{-1} */
      dD_deta_gr[i-1]=Dprime_gr[i-1]*a*a*Hdim;                     /* dD/deta */
      //printf("%d %f %f %f\n", i, aa_gr[i], xx_gr[i], Dprime_gr[i]);
  }

   my_spline(xx_gr, DD_gr, N_GROWTH, 0, 1e31,y2_DD);
   my_spline(xx_gr, gg_gr, N_GROWTH, 0, 1e31,y2_gg);
   my_spline(xx_gr, Dprime_gr, N_GROWTH, 1.0, 1e31,y2_Dprime);
   my_spline(xx_gr, dD_deta_gr, N_GROWTH, 1.0, 1e31,y2_dD_deta);

 /* return dlnD/dlna, really */
  return(deldot_return);
}
void print_matrix(double **a, int n)
{
  int i, j;

  for(i=1; i<=n; i++)
    {
      for(j=1; j<=n; j++)
	printf("%e ", a[i][j]);
      printf("\n");
    }
  printf("\n"); 
}

void file_print_matrix(double **a, int n, char *file)
{
  int i, j;
  FILE *ifp;

  ifp=fopen(file, "w");
  for(i=1; i<=n; i++)
    {
      for(j=1; j<=n; j++)
          fprintf(ifp, "%e ", a[i][j]);
      fprintf(ifp, "\n");
    }
  fclose(ifp);
}
void file_read_matrix(double **a, int n, char *file)
{
  int i, j;
  FILE *ifp;

  ifp=fopen(file, "r");
  for(i=1; i<=n; i++)
    {
      for(j=1; j<=n; j++)
          fscanf(ifp, "%lf ", &a[i][j]);
      fscanf(ifp, "\n");
    }
  fclose(ifp);
}
void print_matrix_nonsquare(double **a, int m, int n)
{
  int i, j;

  for(i=1; i<=m; i++)
    {
      for(j=1; j<=n; j++)
	printf("%e ", a[i][j]);
      printf("\n");
    }
  printf("\n"); 
}
void print_histogram(gsl_histogram *h, char *filename)
{
  size_t bins, i;
  double lower, upper;
  FILE *ifp;

  bins=gsl_histogram_bins(h);

  ifp=fopen(filename, "w");
  for(i=0; i < bins; i++) 
    {
      gsl_histogram_get_range(h, i, &lower, &upper); 

      fprintf(ifp, "%f %f\n",  lower, gsl_histogram_get(h, i));
      fprintf(ifp, "%f %f\n",  upper, gsl_histogram_get(h, i));
      /**printf("%d %f\n", i, gsl_histogram_get(h, i));*/
    }
  fclose(ifp);
}
double inverse_lndet_gsl_LU(double **A, double **Ainv, int n)
{
    /************************************/
    /* same, but use the GSL routines  **/
    /* as a bonus, returns ln(det(A))  **/
    /************************************/
    int s, i, j;
    gsl_matrix * a=gsl_matrix_alloc(n,n), * ainv=gsl_matrix_alloc(n,n);
    gsl_permutation * p = gsl_permutation_alloc (n);
    double lndet;
    
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            gsl_matrix_set (a, i, j, A[i+1][j+1]);
    
    gsl_linalg_LU_decomp (a, p, &s);    // a is turned into LU
    gsl_linalg_LU_invert (a, p, ainv);  // takes in LU
    
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            Ainv[i+1][j+1]=gsl_matrix_get (ainv, i, j);

    lndet = gsl_linalg_LU_lndet (a);  // log(Det(LU))
    
    gsl_permutation_free (p);
    gsl_matrix_free(a) ;
    gsl_matrix_free(ainv);

    return(lndet);
}
void inverse(double **a, double **ai, int n)
/* takes a matrix **a (nxn) and finds its inverse and stores it in **ai */
{
       int k, l;
       double  **b,  **x;
       x=dmatrix(1,n, 1, 1);
       b=dmatrix(1,n,1,1);
       /* save matrices for later testing of results */
       for (l=1;l<=n;l++) {
           for (k=1;k<=n;k++) ai[k][l]=a[k][l];
           for (k=1;k<=1;k++) x[l][k]=b[l][k];
       }
   /* invert matrix */
       my_gaussj(ai,n,x,1);
    /* this is a NR routine which does the job*/             

       free_dmatrix(x,1,n, 1, 1);
       free_dmatrix(b,1,n, 1, 1);
}
void mat_mat_prod(double **a, double **b, double **c, int dim)
     /* multiplies matrix **a and matrix **b and stores it matrix **c */
{     
  double sum;
  int i,j,k;
  for(i=1; i<=dim; i++)
    for(k=1; k<=dim; k++)
      {
	sum=0.0;
	for(j=1; j<=dim; j++)
	  sum+=a[i][j]*b[j][k];
	c[i][k]=sum;
      }
}
int ij_to_alpha(int i, int j)
{
    /*************************************************************************/
    /** convert from redshift-bins i and j to single variable alpha         **/
    /** assumes that i goes from 1 to Z_BINS, and j goes from i to Z_BINS-1 **/
    /*************************************************************************/
    return (Z_BINS*(i-1) - i*(i-1)/2 + j);
}

double Hsq(double z, double omega_m, double eqstate, double wprime)
/*******************************/
/*** in h/Mpc (squared) !!!  ***/
/*******************************/
{

  return H0_hmpcinv*H0_hmpcinv*(omega_m*pow(1+z, 3) + 
                (1-omega_m)*pow(1+z, 3*(1+eqstate+wprime))*exp(-3*wprime*z/(1+z)) );
}
double H_dimless(double z, double omega_m, double eqstate, double wprime)
/*******************************/
/*** in h/Mpc (squared) !!!  ***/
/*******************************/
{

    return sqrt(omega_m*pow(1+z, 3) + 
                (1-omega_m)*pow(1+z, 3*(1+eqstate+wprime))*exp(-3*wprime*z/(1+z)) );
}
double D_comov(double zmax, double omega_m, double eqstate, double wprime)
{
    /*******************************/
    /*** in Mpc/h  ***/
    /*******************************/
    double z, sum, dz;

  dz=0.001;
  sum=0; 
  for(z=0; z<=zmax-VERYSMALL; z+=dz)
      sum+=dz/sqrt(Hsq(z, omega_m, eqstate, wprime));

  return sum;

}
double nz(double z)
{   
    return(0.5*pow(Z_0, -3)*z*z*exp(-z/Z_0));
}    
