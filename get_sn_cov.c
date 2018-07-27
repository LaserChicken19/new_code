#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<omp.h>
#include"recipes/nrutil.h" 
#include"recipes/nr.h"
#include "declarations.h" 

#define BETA_CMB_DIP (370.0/3.0e5)
#define VEC_X        -0.071 // dipole x, y, z
#define VEC_Y        -0.662
#define VEC_Z         0.746

void read_SN_pos(double **SN_z_th_phi, char *filename)
{
    int i;
    int idum;
    double z,delta_m, sigma_m, ra, dec;
    char sdum[256];
    FILE *ifp1;

    N_SN = file_eof_linecount(filename);
    
    ifp1=fopen(filename, "r");
    if (ifp1 == NULL) {
        fprintf(stderr, "Can't open input file in read_SN_pos!\n");
        exit(1);
    }
    
    for(i=1; i<=N_SN; i++)
    {
        fscanf(ifp1, "%d %s %lf %lf %lf %lf %lf\n",
               &idum, sdum, &z,&delta_m,&sigma_m, &ra, &dec);
        SN_z_th_phi[i][1] = z;
        SN_z_th_phi[i][2] = 90.0-dec;
        SN_z_th_phi[i][3] = ra;
        SN_z_th_phi[i][4] = sigma_m;
    }
    fclose(ifp1);
}
void read_pos_noise(int OFFDIAG, double **SN_z_th_phi, double *delta_m, char *filename)
{
    int i,j;
    double z,m,cov,sigma_m, m_th, l_gal, b_gal;

    FILE *ifp1;
    N_SN = file_eof_linecount(filename);

    ifp1=fopen(filename, "r");
    if (ifp1 == NULL) {
        fprintf(stderr, "Can't open input SN file\n");
        exit(1);
    }
    
    for(i=1; i<=N_SN; i++)
    {
        fscanf(ifp1, "%lf %lf %lf %lf %lf %lf\n",
               &z, &m, &sigma_m, &l_gal, &b_gal, &m_th);
        
        SN_z_th_phi[i][1] = z;
        SN_z_th_phi[i][2] = (90.0-b_gal) * M_PI/180.0;
        SN_z_th_phi[i][3] = l_gal* M_PI/180.0;
        SN_z_th_phi[i][4] = sigma_m;
        delta_m[i]        = m - m_th;

        /***********************************************************************/
        /** OPTIONALLY, MAKE CORRECTION SUGGESTED BY CASTRO - 2nd term only:  **/
        /** It leads to very small (but noticeable - thickness of thick line) **/
        /** differences in likelihoods of A and v_extra - Dec 2015            **/
        /***********************************************************************/
        /*
   double xsn, ysn, zsn, corr;
        xsn = sin(SN_z_th_phi[i][2])*cos(SN_z_th_phi[i][3]);
        ysn = sin(SN_z_th_phi[i][2])*sin(SN_z_th_phi[i][3]);
        zsn = cos(SN_z_th_phi[i][2]);
        corr = 5*log10(1 - BETA_CMB_DIP * (xsn*VEC_X + ysn*VEC_Y + zsn*VEC_Z));
        delta_m[i] -= corr;
        */
    }
    fclose(ifp1);

    //if (OFFDIAG != -1)
   // {
     //   ifp1=fopen(file_noise_cov, "r");
     //   for(i=1; i<=N_SN; i++)
     //   {
     //       for(j=1; j<=N_SN; j++)
     //       {
     //           fscanf(ifp1, "%lf ", &cov);
                //if ((i == 208 || i == 207) && (j == 208)) printf("Cov=%f\n", cov);
     //           Noise_Cov[i][j] = cov;
     //       }
     //       fscanf(ifp1, "\n");
     //   }
     //   fclose(ifp1);
   // }
   // else // if using N_ii = sigma^2        
   // {
     //   for(i=1; i<=N_SN; i++)
      //      for(j=1; j<=N_SN; j++)
       //     {
        //        if (i == j)
        //        {
         //           z = SN_z_th_phi[i][1];
                    //prefac = 5/log(10.0) * (1 - (1+z)/H_dimless(z, 0.3, -1.0, 0.0)/r_dimless_tab(z));
                    // no need to add 300km/s any more
                    //Noise_Cov[i][j] = pow(SN_z_th_phi[i][4], 2) + pow(prefac * 300.0/3.0e5, 2);
           //         Noise_Cov[i][j] = pow(SN_z_th_phi[i][4], 2);
            //    }
          //      else        Noise_Cov[i][j] = 0;
           // }
   // }
}
void order_SN_increasing_z(int nsn, double **SN_z_th_phi, double **Noise_Cov, double *delta_m)
{
    /**********************************************************************************/
    /** orders all of the elements of SN_z_th_phi[nsn][4] in order of increasing z, **/
    /** which is the [*][1] element. Returns the reordered multi-array SN_z_th_phi   **/
    /***********************************************************************************/
    int i, j;
    double *zarr, **sn_mat_copy, *delta_m_copy, **Noise_Cov_copy;
    unsigned long *indx;
    FILE *ifp;
    zarr = dvector(1, nsn);
    indx = lvector(1, nsn);
    delta_m_copy  = dvector(1, nsn);
    sn_mat_copy   = dmatrix(1, nsn, 1, 4);
    Noise_Cov_copy= dmatrix(1, nsn, 1, nsn);

    for (i=1; i<=nsn; i++)
    {
        zarr[i] = SN_z_th_phi[i][1];

        delta_m_copy[i] = delta_m[i];
        for(j=1; j<=4; j++)
            sn_mat_copy[i][j] = SN_z_th_phi[i][j];
        for(j=1; j<=nsn; j++)
            Noise_Cov_copy[i][j] = Noise_Cov[i][j];    
    }
    
    my_indexx(nsn, zarr, indx);

    for (i=1; i<=nsn; i++)
    {
        delta_m[i] = delta_m_copy[indx[i]];

        for(j=1; j<=4; j++)
            SN_z_th_phi[i][j] = sn_mat_copy[indx[i]][j];

        for(j=1; j<=nsn; j++)
            Noise_Cov[i][j]   = Noise_Cov_copy[indx[i]][indx[j]];
    }
    ifp=fopen("ordered_z_SN.dat", "w");
    for(i=1; i<=nsn; i++)
        fprintf(ifp, "%d %f\n", i, SN_z_th_phi[i][1]);
    fclose(ifp);
}
void order_whole_SN_file_in_z(char *file_in, char *file_out)
{
  
    int i,nsn;
    double z,m,sigma_m, m_th, l_gal, b_gal;
    double *zarr, **arr;
    unsigned long *indx;
    FILE *ifp;
    nsn = file_eof_linecount(file_in);

    arr = dmatrix(1, nsn, 1, 6);
    zarr= dvector(1, nsn );
    indx = lvector(1, nsn);

    ifp=fopen(file_in, "r");
    if (ifp == NULL) {
        fprintf(stderr, "Can't open input file in  order_whole_SN_file_in_z\n");
        exit(1);
    }

    for(i=1; i<=nsn; i++)
    {
        fscanf(ifp, "%lf %lf %lf %lf %lf %lf\n",
               &z, &m, &sigma_m, &l_gal, &b_gal, &m_th);

        zarr[i]  = z;
        arr[i][1] = z;
        arr[i][2] = m;
        arr[i][3] = sigma_m;
        arr[i][4] = l_gal;
        arr[i][5] = b_gal;
        arr[i][6] = m_th;

    }
    fclose(ifp);

    my_indexx(nsn, zarr, indx);

    ifp=fopen(file_out, "w");
    for(i=1; i<=nsn; i++)
        fprintf(ifp, "%f %f %f %f %f %f\n",
                arr[indx[i]][1], arr[indx[i]][2], arr[indx[i]][3],
                arr[indx[i]][4] ,arr[indx[i]][5], arr[indx[i]][6]);
    fclose(ifp);
}
    



void calculate_Cov_vel_of_SN(int i_1,int j_1, int ij_size,  double **SN_z_th_phi, double **Signal_SN,
                             double omega_m, double w0, double wa)
{

    int  i;
#pragma omp parallel for default(none) shared(i_1,j_1,ij_size,  Signal_SN, SN_z_th_phi, omega_m, w0, wa) schedule(dynamic)  
    
    for(i=i_1; i<=i_1+ij_size-1; i++)
    {
        //int th_id = omp_get_thread_num();
        //printf("i=%d from thread  %d, redshift=%f\n", i, th_id, SN_z_th_phi[i][1]);
        
        int j;
        double cosalpha, dum;    

        for(j=j_1+i-i_1; j<=j_1+ij_size-1 ; j++) // upper triangular
        {
            if (i == j) cosalpha=1.0;
            else 
                cosalpha = getang_spher_cos(SN_z_th_phi[i][2], SN_z_th_phi[i][3],
                                            SN_z_th_phi[j][2], SN_z_th_phi[j][3]);

            dum = Cij_theta_gsl_int(i, j, SN_z_th_phi[i][1], SN_z_th_phi[j][1], 
                                    cosalpha, omega_m, w0, wa);
            //printf("calculating i=%d, j=%d\n",i-i_1+1,j-j_1+1);
            Signal_SN[i-i_1+1][j-j_1+1] = dum;
            Signal_SN[j-j_1+1][i-i_1+1] = dum;

        }
    }
}
void calculate_Signal_given_z_theta_arr(int ncosalpha, int nz, double *cosalpha_arr, double *z_arr, double ***Signal_tensor,
                                        double omega_m, double w0, double wa)
{

    int  ncos;
#pragma omp parallel for default(none) shared(ncosalpha, nz, Signal_tensor, cosalpha_arr, z_arr, omega_m, w0, wa) schedule(dynamic)  
    
    for(ncos=1; ncos<=ncosalpha; ncos++)
    {
        int th_id = omp_get_thread_num();
        ///printf("i=%d from thread  %d cosalpha=%f\n", ncos, th_id, cosalpha_arr[ncos]);

        double cosalpha = cosalpha_arr[ncos];
        int i;
        for(i=1; i<=nz; i++)
        {
            int j;
            double dum;               
            for(j=i; j<= nz; j++) // upper triangular
            {
                dum = Cij_theta_gsl_int(i, j, z_arr[i], z_arr[j], cosalpha, omega_m, w0, wa);
                
                Signal_tensor[ncos][i][j] = dum;
                Signal_tensor[ncos][j][i] = dum;
                
            }
        }
    }
}   

double getang_spher_cos(double theta1, double phi1, double theta2, double phi2)
     /* compute the angle between the first and second point */
{
  double x1, y1, z1, x2, y2, z2;

  x1=sin(theta1)*cos(phi1);
  y1=sin(theta1)*sin(phi1);
  z1=cos(theta1);

  x2=sin(theta2)*cos(phi2);
  y2=sin(theta2)*sin(phi2);
  z2=cos(theta2);

  return(x1*x2 + y1*y2 + z1*z2);
}
double trace(double **a, int dim)
/* takes trace of a matrix **/
{     
  double sum;
  int i;

  sum=0.0;
  for(i=1; i<=dim; i++)
    sum+=a[i][i];

  return(sum);
}
double Fisher(double **Signal, double **Noise)
{
    /*******************************************/
    /*** F = 1/2 Trace[(C^{-1}C,A)^2] where  **/
    /*** C = A*S+N with fiducial A=1         **/
    /*** C,A = S                             **/
    /******************************************/
    int i, j;
    double **Cinv, **dum1, **dum3, **Ctot;

    Cinv = dmatrix(1, N_SN, 1, N_SN);
    Ctot = dmatrix(1, N_SN, 1, N_SN);
    dum1 = dmatrix(1, N_SN, 1, N_SN);
    dum3  = dmatrix(1, N_SN, 1, N_SN);

    for(i=1; i<=N_SN; i++)
        for(j=1; j<=N_SN; j++)
            Ctot[i][j] = Signal[i][j] + Noise[i][j];

    inverse(Ctot, Cinv, N_SN);

    //  C^{-1} C_A
    mat_mat_prod(Cinv, Signal, dum1, N_SN);

    // [C^{-1} C_A]^2
    mat_mat_prod(dum1, dum1, dum3, N_SN);

    ///print_matrix(Signal, N_SN);

    free_dmatrix(Cinv, 1, N_SN, 1, N_SN);
    free_dmatrix(Ctot, 1, N_SN, 1, N_SN);
    free_dmatrix(dum1, 1, N_SN, 1, N_SN);

    return(0.5*trace(dum3, N_SN));
}

void check_condition_number_matrix(int n, double **matrix)
{
    // finds condition number of matrix M, but ignores to take abs value of eigenvalues
    int nrot, i , j;
    double **eigenvec, *eigenval, **mat;
    
    mat        = dmatrix(1, n, 1, n);
    eigenvec   = dmatrix(1, n, 1, n);
    eigenval   = dvector(1, n);

    // save the original first
    for (i=1; i<=n; i++)
        for (j=1; j<=n; j++)
            mat[i][j]=matrix[i][j];
    
    my_jacobi(mat, n, eigenval, eigenvec, &nrot);    
    my_eigsrt(eigenval, eigenvec, n);

    //printf("Eigenvalues are\n");
    // for (i=1; i<=n; i++)
    //     printf("%d %e\n", i, eigenval[i]);
    //printf("The first eigenvector is:\n");
    //for (i=1; i<=n; i++)
    //    printf("%e ", eigenvec[i][1]);
    //printf("\n");

    printf("The condition number is %e\n", eigenval[1]/eigenval[n]);
}
double log_det_matrix(int n, double **matrix)
{
    /*************************************************************************/
    /** calculate the log of determinant as a sum over the log(eigenvalues) **/
    /*************************************************************************/
    int nrot, i , j;
    double **eigenvec, *eigenval, **mat, log_det=0;
    
    mat        = dmatrix(1, n, 1, n);
    eigenvec   = dmatrix(1, n, 1, n);
    eigenval   = dvector(1, n);

    // save the original first
    for (i=1; i<=n; i++)
        for (j=1; j<=n; j++)
            mat[i][j]=matrix[i][j];
    
    my_jacobi(mat, n, eigenval, eigenvec, &nrot);    
    my_eigsrt(eigenval, eigenvec, n);

    //printf("Eigenvalues are\n");
    for (i=1; i<=n; i++)
        printf("%d %e\n", i, eigenval[i]);

    log_det = 0;
    for (i=1; i<=n; i++)
        log_det +=  log(eigenval[i]);
    return(log_det);
    
}

double chisq_cov(int n, double A, double *mu, double **Signal, double **Noise)
{
    int i=0, j=0;
    double **Cov,  **Covinv;
    double first_term=0, second_term=0;
    
    Cov    = dmatrix(1, n, 1, n);
    Covinv = dmatrix(1, n, 1, n);

    
    for(i=1; i<=n; i++)
        for(j=1; j<=n; j++)
            Cov[i][j] = A*Signal[i][j] + Noise[i][j];

    inverse(Cov, Covinv, n);

    for (j=1; j<=n; j++)
        for (i=1; i<=n; i++)
            first_term += mu[i]*Covinv[i][j]*mu[j];

    second_term = log_det_matrix(n, Cov);
    //printf("first term=%f, second term(log_det)=%f\n", first_term, second_term);
    
    return(first_term+second_term);
}
double chisq_cov_input_Cinv(int iter, int n, double *mu, double ***Cov_3D, double *logdetCov_arr, double ***Covinv_3D)
{
    /////////////////////////////////////////////////////////////////////////
    // same as chisq_cov, but to save time you pre-computed Cinv using openmp
    /////////////////////////////////////////////////////////////////////////
    int j, k;
    double first_term=0, second_term=0;
    double **Cov, **Covinv, *dm;
    
    // first, map the 3D arrays into the desired 2Ds
    Cov    = dmatrix(1, n, 1, n);
    Covinv = dmatrix(1, n, 1, n);
    dm     = dvector(1, n);

    // if so desired, covert mag into an array of \tilde(\mu) which account for
    // fact that the likelihood is really gaussian in mu, not mag - see Dragan's note
    // from April 2016
    if (NG_LIKE_CORRECTION == 1)
        for(j=1; j<=n; j++)
            dm[j] = 5.0 / log(10.0) * (pow(10.0, 0.2*mu[j]) - 1); // mag(rhs) -> mu(lhs)
    else
        for(j=1; j<=n; j++)
            dm[j] = mu[j]; // leave magnitudes 

    for(j=1; j<=n; j++)
        for(k=1; k<=n; k++)
        {
            Cov   [j][k] = Cov_3D   [iter][j][k];
            Covinv[j][k] = Covinv_3D[iter][j][k];
        }
 
    
    for (j=1; j<=n; j++)
        for (k=1; k<=n; k++)
            first_term += dm[j]*Covinv[j][k]*dm[k];
    
    /******************************************************************/
    /** if dm=magnitudes, marg over scriptM is just the usual formula */
    /******************************************************************/
    if ((MARG_OVER_SCRIPTM == 1) && (NG_LIKE_CORRECTION == 0))
    {
        double X10=0, X11=0;
        
        for (j=1; j<=n; j++)
            for (k=1; k<=n; k++)
            {
                X10 += 1.0*Covinv[j][k]*dm[k];
                X11 += 1.0*Covinv[j][k]*1.0;
            }
        // first_term = X00 + log(X11/2pi) - X10^2/X11;
        first_term += log10(X11/(2*M_PI)) - X10*X10/X11;
    }

    /******************************************************************/
    /** if dm = magnification, need to numerically marg over scriptM **/
    /** since a flat prior in scriptM is a different prior in mu    ***/
    /******************************************************************/
    if ((MARG_OVER_SCRIPTM == 1) && (NG_LIKE_CORRECTION == 1))
        first_term = chisq_in_mu_marg_over_scriptM(n, dm, Covinv, -0.2, 0.2, 0.02);

    second_term = logdetCov_arr[iter] ;

    free_dmatrix(Cov,    1, n, 1, n);
    free_dmatrix(Covinv, 1, n, 1, n);
    free_dvector(dm, 1, n);
    
    return(first_term+second_term);
}
double chisq_in_mu_marg_over_scriptM(int n, double* v, double **Cinv, double Mmin, double Mmax, double dM)
{
    // input chisq = v^T Cinv v, and min, max and step in scriptM
    // assumes v = vector of magnifications mu (for mags this is analytic and not using this fun)
    // output -2*ln(likelihood_from_chisq_marg_over_scriptM)

    int  j, k;
    double chisq=0, M, L=0, fac;
    double *v1;

    v1 = dvector(1, n);

    
    for(M=Mmin; M<=Mmax+1.0e-10; M+=dM)
    {
        // get offsets
/*
        if (NG_LIKE_CORRECTION == 0) // Gaussian in magnitudes
            for (j=1; j<=n; j++)
                v1[j] = v[j] + M;
*/      
//        if (NG_LIKE_CORRECTION == 1) // Gaussian in magnification - NOW DEFAULT HERE
//        {
            fac = pow(10.0, 0.2*M);
            for (j=1; j<=n; j++)
                v1[j] = fac * (v[j]+5.0/log(10.0)) - 5.0/log(10.0);
//        }
        
        chisq = 0;
        for (j=1; j<=n; j++)
            for (k=1; k<=n; k++)
                chisq += v1[j]*Cinv[j][k]*v1[k];

        L += exp(-0.5*chisq) * dM;
    }

    L = L /(Mmax-Mmin);

    free_dvector(v1, 1, n);
    
    return(-2.0*log(L));
}        
               
void find_hist_68_95_99_CL(gsl_histogram *h_orig, double *peak_min_max)
{
    /******************************************************************/
    /** In a 1-D histogram, find max-L model and the confid contours***/
    /** peak_max_min[1]   is the highest-likelihood value           ***/
    /** peak_max_min[2,3] are the lowest and highest at 68.3% CL    ***/
    /** peak_max_min[4,5] are the lowest and highest at 95.4% CL    ***/
    /** peak_max_min[6,7] are the lowest and highest at 99.7% CL    ***/
    /*******************************************************************/
    int j=0;
    double total, sum=0;    
    double lower, upper, v, val, frac;
    gsl_histogram *h;

    h = gsl_histogram_clone(h_orig);
    
    // get the maximum-likelihood x-value 
    j = gsl_histogram_max_bin (h);                // bin corresponding to highest value
    gsl_histogram_get_range(h,j,&lower,&upper);   // v(k) (bin boundary) values at that value
    v = 0.5*(lower+upper);
    peak_min_max[1] = v;

    // initialize these to ridiculously low and high values
    peak_min_max[2] = peak_min_max[4] = peak_min_max[6] = +1.0e20;
    peak_min_max[3] = peak_min_max[5] = peak_min_max[7] = -1.0e20;

    // find total sum in histogram
    total = gsl_histogram_sum (h);

    // for each confidence level, find lower and upper bound in x
    frac = 0.683;    
    while (sum < frac*total)
    {
        val = gsl_histogram_max_val (h);              // highest current value in histogram
        sum += val;            
        j = gsl_histogram_max_bin (h);                // bin corresponding to highest value
        gsl_histogram_get_range(h,j,&lower,&upper);   // v (bin boundary) values at that value
        v = 0.5*(lower+upper);
        gsl_histogram_accumulate(h, v, -val);       // zero out that histogram bin

        if (v < peak_min_max[2])  peak_min_max[2]=v;
        if (v > peak_min_max[3])  peak_min_max[3]=v;
    }
    printf("Maximum value = %f in hist bin %d\n", peak_min_max[1], j);
    //printf("Low  68.3 bound = %f\n", peak_min_max[2]);
    //printf("High 68.3 bound = %f\n", peak_min_max[3]);

    // inherit the previous boundaries
    // (crucial of you already hit a wall at one side!)
    peak_min_max[4] = peak_min_max[2];
    peak_min_max[5] = peak_min_max[3];
    
    frac = 0.954;    
    while (sum < frac*total)
    {
        val = gsl_histogram_max_val (h);              // highest current value in histogram
        sum += val;            
        j = gsl_histogram_max_bin (h);                // bin corresponding to highest value
        gsl_histogram_get_range(h,j,&lower,&upper);   // v(k) (bin boundary) values at that value
        v = 0.5*(lower+upper);
        gsl_histogram_accumulate(h, v, -val);       // zero out that histogram bin

        if (v < peak_min_max[4])  peak_min_max[4]=v;
        if (v > peak_min_max[5])  peak_min_max[5]=v;
    }
    //printf("Low  95.3 bound = %f\n", peak_min_max[4]);
    //printf("High 95.3 bound = %f\n", peak_min_max[5]);

    // inherit the previous boundaries
    // (crucial of you already hit a wall at one side!)
    peak_min_max[6] = peak_min_max[4];
    peak_min_max[7] = peak_min_max[5];

    frac = 0.997;    
    while (sum < frac*total)
    {
        val = gsl_histogram_max_val (h);              // highest current value in histogram
        sum += val;            
        j = gsl_histogram_max_bin (h);                // bin corresponding to highest value
        gsl_histogram_get_range(h,j,&lower,&upper);   // v(k) (bin boundary) values at that value
        v = 0.5*(lower+upper);
        gsl_histogram_accumulate(h, v, -val);       // zero out that histogram bin

        if (v < peak_min_max[6])  peak_min_max[6]=v;
        if (v > peak_min_max[7])  peak_min_max[7]=v;
    }
    //printf("Low  99.7 bound = %f\n", peak_min_max[6]);
    //printf("High 99.7 bound = %f\n", peak_min_max[7]);

}



double calculate_norm_like_A(int nsn,
                             double Amin, double Amax, double dA,
                             double sig_int_sq_min, double sig_int_sq_max, double dsig_int_sq,
                             double *deltam, double **Signal, double **Noise, char *filename)
{
    // returns the p-value too!
    
    int i, N_A, N_SIG_INT_SQ, NTOT;
    double *A_arr, *L_arr, max_like, min_chi=1.0e10, A_gt_1_like=0, sum=0;
    //char filename[256];
    FILE *ifp;

    //sprintf(filename, "%s%d%s", file_root, nsn, ".dat");
    printf("filename %s\n", filename);

    N_A          = round((Amax-Amin)/dA) + 1;
    N_SIG_INT_SQ = round((sig_int_sq_max-sig_int_sq_min)/dsig_int_sq) + 1;
    NTOT         = N_A * N_SIG_INT_SQ;
        printf("N_A, N_SOG_INT_SQ, NTOT are %d %d %d\n", N_A, N_SIG_INT_SQ, NTOT);
    
    L_arr = dvector(1, N_A);
    A_arr = dvector(1, N_A);

    /******************************************************************/
    /*** new - precompute the inverse of the A*S+N using openmp     ***/
    /*** note, logdet(Cov) is a HUGE time hog; parallelize that, too **/
    /******************************************************************/
    double ***Cov_3D, *logdetCov_arr;
    double ***Covinv_3D;
    Cov_3D    = d3tensor(1, NTOT, 1, nsn, 1, nsn);
    Covinv_3D = d3tensor(1, NTOT, 1, nsn, 1, nsn);
    logdetCov_arr = dvector(1, NTOT);

#pragma omp parallel for default(none) shared(nsn, Amin, dA, sig_int_sq_min, dsig_int_sq, N_A, N_SIG_INT_SQ, NTOT, Signal, Noise, Cov_3D, logdetCov_arr, Covinv_3D) schedule(dynamic) 
    
    for(i=1; i<=NTOT; i++)
    {
        int i1, i2;
        
        i_to_i1i2(i, &i1, &i2, N_A, N_SIG_INT_SQ);

        double A = Amin + (i1-1)*dA;
        double sig_int_sq = sig_int_sq_min + (i2-1)*dsig_int_sq;
        
        //int th_id = omp_get_thread_num();
        //printf("i=%d, A=%f from thread  %d\n", i, A, th_id);
        
        int j, k;
        double **Cov, **Covinv;
        Cov    = dmatrix(1, nsn, 1, nsn);
        Covinv = dmatrix(1, nsn, 1, nsn);
    
        for(j=1; j<=nsn; j++)
            for(k=1; k<=nsn; k++)
            {
                Cov[j][k] = A*Signal[j][k] + Noise[j][k];
                if (j == k)
                    Cov[j][k] += sig_int_sq;
            }

        /*******************************************/
        /* calculate ln(Det) and inverse in one go */
        /*******************************************/

        logdetCov_arr[i] = inverse_lndet_gsl_LU(Cov, Covinv, nsn);
        //inverse(Cov, Covinv, nsn);
        //logdetCov_arr[i] = log_det_matrix(nsn, Cov);
        
        for(j=1; j<=nsn; j++)
            for(k=1; k<=nsn; k++)
            {
                Cov_3D   [i][j][k] = Cov   [j][k];
                Covinv_3D[i][j][k] = Covinv[j][k];
            }
    
        free_dmatrix(Cov,    1, nsn, 1, nsn);
        free_dmatrix(Covinv, 1, nsn, 1, nsn);
    }
    printf("Done with Cinv calculations\n");

    /***********************************************************************/
    /*** now do the usual and calculate likelihood given Cinv for each A ***/
    /***********************************************************************/
    int i1;
#pragma omp parallel for default(none) shared(nsn, Amin, dA, A_arr, N_A, L_arr, deltam, N_SIG_INT_SQ, Cov_3D, logdetCov_arr, Covinv_3D, min_chi, max_like, sum, A_gt_1_like) schedule(dynamic)

    for(i1=1; i1<=N_A; i1++)    
    {
        double A = Amin + (i1-1)*dA;
        
        A_arr[i1] = A;
        L_arr[i1] = 0; //initialize

        if (fmod(A_arr[i1], 0.5) < 1.0e-5)
            printf("Working on A=%f\n", A_arr[i1]);
        
        // OLDEST, NO PARALLEL: chi = chisq_cov(nsn, A, deltam, Signal, Noise);
        // OLD: NO MARG OVER SIG_INT: chi = chisq_cov_input_Cinv(i1, nsn, deltam, Cov_3D, logdetCov_arr, Covinv_3D);
        int i2, ii;
        double chii;
        for(i2=1; i2<=N_SIG_INT_SQ; i2++)
        {
            ii = (i1-1) * N_SIG_INT_SQ + i2;
            chii = chisq_cov_input_Cinv(ii, nsn, deltam, Cov_3D, logdetCov_arr, Covinv_3D);
            
            // temporary, check how L(A=1, sig_int^sq) depends on the latter
            // printf("%f %e\n", sig_int_sq_min + (i2-1)*dsig_int_sq, exp(-0.5*chi)/3.425020e+91); 
            L_arr[i1] += exp(-0.5*chii);            
        }

        // overall chi, marg over sigma_int^{sq}
        chii = -2.0*log(L_arr[i1]);
        
        if (chii < min_chi)
        {
            min_chi = chii;
            max_like = L_arr[i1];
        }

        // for p-value
        sum += L_arr[i1];
        if (A > 1.0) A_gt_1_like += L_arr[i1];
    }
    printf("Done with A likelihood, sum=%e\n", sum);
    
    /***********************************************************/
    /**** create the histogram of A to get 68%, 95% CL ranges **/
    /***********************************************************/
    gsl_histogram *Ahist;
    double *peak_min_max;
    peak_min_max = dvector(1, 7);
    
    Ahist = gsl_histogram_alloc (N_A);  // each round value is at center of bin
    gsl_histogram_set_ranges_uniform (Ahist, Amin-0.5*dA, Amax+0.5*dA);

    for(i1=1; i1<=N_A; i1++)
        gsl_histogram_accumulate(Ahist, A_arr[i1], L_arr[i1]);
    //print_histogram(Ahist, "OUTPUT_PAPER2/test_hist.dat");
    
    find_hist_68_95_99_CL(Ahist, peak_min_max);
    printf("************** likelihood range analysis of A *************\n");
    printf("Highest-L A: %f\n",   peak_min_max[1]);
    printf("68perc CL range: [%f, %f]\n", peak_min_max[2], peak_min_max[3]);
    printf("95perc CL range: [%f, %f]\n", peak_min_max[4], peak_min_max[5]);

    /***************************************/
    /**** report a couple of Delta chi^2  **/
    /***************************************/
    i1=1;
    printf("-2*Delta ln(L) for bestfit vs A=%f: %f\n", A_arr[i1], -2*log(L_arr[i1]/max_like));
    i1 = round(1.0/dA+1);
    printf("-2*Delta ln(L) for bestfit vs A=%f: %f\n", A_arr[i1], -2*log(L_arr[i1]/max_like));
    
    /************************/
    /**** write to a file  **/
    /************************/
    ifp=fopen(filename, "w");
    for(i1=1; i1<=N_A; i1++)
        fprintf(ifp, "%f %f\n", A_arr[i1], L_arr[i1]/max_like);
    fclose(ifp);

    return(A_gt_1_like/sum);
}

void i_to_i1i2(int i, int *i1, int *i2, int i1max, int i2max)
{

    ////////////////////////////////////////////////////////////////////
    // given composite integer i, returns i1 and i2
    // all indices start at 1
    // relation is i = (i1-1) * i2max + (i2-1) + 1, where i1 = [1,i1max], i2 = [1,i2max] 
    ////////////////////////////////////////////////////////////////////

    div_t output;

    output = div(i-1, i2max);

    *i2 = output.rem + 1; 
    *i1 = output.quot + 1;
}
