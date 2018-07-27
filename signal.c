#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include<omp.h>
#include"declarations.h"
#include"recipes/nrutil.h"
#include"recipes/nr.h"




int main()
{
    int  OFFDIAG, N_START, N_END, N_USED;
    int i, j, jj;
    double w0, wa, omega_m, A, h, A_k005;
    double  k, n, z;
    //double camb_z=0.1;
    FILE *ifp;
    gsl_set_error_handler_off();
    PRINT_FLAG = 0;
    /**************************************/
    /** Planck parameters - August 2016! **/
    /**************************************/
    omega_m=0.3; 

    n=0.965;  

    omhh=0.140;   
    obhh=0.0223; 

    h=sqrt(omhh/omega_m);
    A_k005 = 2.0e-9; 

    A    = A_k005;

    // only for my own P(k) funcs in functions.c
    A_k0002 = A_k005 / pow(25.0, n-1);

    H0_mpcinv  = h/2997.9;   /* in Mpc-1; need this for ps_linear */
    /********************************/
    /** initialize power spectrum ***/
    /********************************/
    w0=-1.0;
    wa=0.0;
    z=0.0;  
    printf("check point 0\n");
    //double s_1=D_tab(z);
    // printf("sigma8 %0.9e\n",sigma8(0.025));
   // printf("sigma8 %0.9e\n",sigma8(0.075));
   // printf("sigma8 %0.9e\n",sigma8(0.125));
   // printf("sigma8 %0.9e\n",sigma8(0.175));
    /*********************/
    /** Spline the r(z) **/
    /*********************/    
    spline_r(omega_m, w0, wa);
    printf("check point 1 passed\n");
    /************************************************/
    /** Spline the T(k) evaluated using CAMB       **/
    /** do_nonlin=1 HERE, and dn/dlnk=0            **/
    /************************************************/
    //z=0.0;
    run_camb_get_Tk_friendly_format(1, omega_m, omhh, obhh, n, 0.0, A, w0,z);
    printf("check point 2 passed\n");
    spline_Tk_from_CAMB(omega_m, w0, wa, A, n);
    printf("check point 3 passed\n");    
    spline_Pk_from_CAMB   (0, omega_m, w0, wa, A, n);
    printf("check point 4 passed\n");    
    spline_Pk_from_CAMB_NR(0, omega_m, w0, wa, A, n);    
    printf("check point 5 passed\n");
    /************************************************/
    /** Spline the P(k) evaluated at z=0           **/
    /** then later rescale it by the growth factor **/
    /************************************************/

    NSPLINE=0;
    for(k=K_START_INTERP; k<=K_END_INTERP; k+=DLNK_PKSPLINE*k) NSPLINE++;    
   // printf("D(z=0)=%f, D(z=1)=%f, D(z=2)=%f\n", D_tab(0.), D_tab(1.0), D_tab(2.0));
    spline_Pk(z, omega_m, w0, wa, A_k0002, n);

    printf("D(z=0)=%f, D(z=1)=%f, D(z=2)=%f\n", D_tab(0.), D_tab(1.0), D_tab(2.0));  
    printf("the raw sigma8 = %f\n", sigma8(0));
    printf("growth_k0_new  = %f\n", growth_k0_new);
   // printf("f= %f\n",dlnD_dlna(z, omega_m, sqrt(omhh/omega_m), w0, wa));
    /*******************************************************************/
    /**   get #SN and define the SN data and covariance mat arrays   ***/
    /*******************************************************************/

    double **all_SN_z_th_phi, **all_Noise_SN, *all_delta_m;
    double **SN_z_th_phi, **Signal_SN,  **Noise_SN, *delta_m;
    char SN_filename[256];
     // Cov_filename[256];

    //sprintf(SN_filename , "/Users/huterer/research/SUPERNOVAE/SN_ANGPS/SN_DATA/PANSTARRS/SN_data.txt");
   // sprintf(Cov_filename, "/Users/huterer/research/SUPERNOVAE/SN_ANGPS/SN_DATA/PANSTARRS/SN_Cov.txt");
    sprintf(SN_filename , "big_data_clean.txt");
   // sprintf(SN_filename, "bigbigbig.txt");
    //sprintf(Cov_filename, "SN_Cov.txt");
    /************************/
    /*** fiducial choices ***/
    /************************/
    printf("check point 6 passed\n");
    OFFDIAG=1;     /*full off-diag Noise */
    N_START = 1;   /* hardcoded before already */
    N_END=file_eof_linecount(SN_filename); /* How many SN would you actually like to USE? DEFAULT=208, but choice doesntmatter */
    N_END=1500;
    N_USED = N_END-N_START+1;    
    printf("Total number of lines is %d\n",N_END);
    WINDOW = 1;    /* Mean subtracted using int over all-sky (DEFAULT) */

    /*********************************************************************/
    /*** Tabulate the jl'(x) function using both GSL and NR             **/
    /** latter gives the same results, and is required for parallel run **/
    /*********************************************************************/
    LMAX_CIJ_THETA=200; // needs 200-300 for higher-z (z>0.05 or so) convergence even for diag elements - Oct 2016

    spline_jlprime_range_ell   (LMAX_CIJ_THETA);
    printf("check point 7 passed\n"); // can't use for parallel
    spline_jlprime_range_ell_NR(LMAX_CIJ_THETA); // need for parallel runs
    printf("check point 8 passed\n");
    /***********************************************************/
    /***   select datafiles and define some arrays/matrices  ***/
    /***********************************************************/
    N_SN = file_eof_linecount(SN_filename);
    printf("I see %d SN in that file\n", N_SN);
   
    all_SN_z_th_phi = dmatrix(1, N_SN, 1, 4);
    //all_Noise_SN  = dmatrix(1, N_SN, 1, N_SN);
    all_delta_m   = dvector(1, N_SN);
    
    SN_z_th_phi = dmatrix(1, N_USED, 1, 4);
    //Noise_SN  = dmatrix(1, N_USED, 1, N_USED);
    delta_m   = dvector(1, N_USED);

    /***********************************************************/
    /** read the SN from file, possibly with noise covariance **/
    /***********************************************************/

    read_pos_noise(OFFDIAG, all_SN_z_th_phi, all_delta_m, SN_filename);
    printf("check point 9 passed\n");
    /************************************************************************/
    /**                   order the SN in increasing redshift              **/
    /** ... unless if you want SN (ordered) plus gals (ordered anew) case! **/
    /************************************************************************/

    //order_SN_increasing_z(N_SN, SN_z_th_phi, Noise_SN, delta_m);  
    //order_SN_increasing_z(N_SN, all_SN_z_th_phi, all_Noise_SN, all_delta_m);
    //for(i=N_START;i<=N_END;i++){
    //    for(j=1;j<=4;j++)
    //        SN_z_th_phi[i][j]=all_SN_z_th_phi[i][j];
    //}

    //for(i=N_START; i<=N_END; i++)
    //{
    //    for(j=1; j<=4; j++)
    //        SN_z_th_phi[i-N_START+1][j] = all_SN_z_th_phi[i][j];   

    //    for(j=N_START; j<=N_END; j++)
    //        Noise_SN[i-N_START+1][j-N_START+1] = all_Noise_SN[i][j];
    //    delta_m[i-N_START+1]  = all_delta_m[i];   
    //}    

    /**********************************************/
    /** Evaluate the signal matrix of PanStarrs ***/
    /**********************************************/
    Signal_SN = dmatrix(1, N_USED, 1, N_USED);        
    COUNT_IN=0; COUNT_BELOW=0;  COUNT_ABOVE=0;
    double t=omp_get_wtime();
    calculate_Cov_vel_of_SN(1,1501,1500,all_SN_z_th_phi, Signal_SN, omega_m, w0, wa);
    double t1=omp_get_wtime();
    printf("*************************************\n");
    printf("time for covariance  evaluation only=%lf\n", t1-t);    
    printf("*************************************\n");
    printf("%e %e\n", Signal_SN[1][1], Signal_SN[N_USED][N_USED]);
   // printf("%e %e\n", Noise_SN[1][1], Noise_SN[N_USED][N_USED]);
    ifp=fopen("more_error_1500_12.dat", "w");
    for(i=1;i<=N_USED;i++){
        for(j=1;j<=N_USED;j++)
            fprintf(ifp,"%0.9e ",Signal_SN[i][j]);
        fprintf(ifp,"\n");
    }
    fclose(ifp);
   // ifp=fopen("noiseMatrix.dat", "w");
   // for(i=1;i<=N_USED;i++){
   //     for(j=1;j<=N_USED;j++)
   //         fprintf(ifp,"%0.9e ",Noise_SN[i][j]);
   //     fprintf(ifp,"\n");
   // }
   // fclose(ifp);
    printf("complete!\n"); 
    exit(0);
} 
