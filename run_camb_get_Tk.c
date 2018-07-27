/************************************************************************/
/** Run CAMB, extract the transfer function, and write it into a file  **/
/** Spline the transfer function in the main code, use it in TF_general**/
/************************************************************************/
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "recipes/nrutil.h" 
#include"recipes/nr.h"
#include"declarations.h"

void run_camb_get_Tk_friendly_format(int do_nonlinear, double omega_m, double omhh_local, double obhh_local, double n, double dn_dlnk, double A_k005, double w0,  double r_s)
{
    int  i;
    double h, H0, ochh;
    double k, Tk_CDM, Tk_bar, Tk_phot, Tk_nu_massless, Tk_nu_massive, Tk_tot, Tk_norm;
    FILE *ifp , *ifp1;

    h=sqrt(omhh_local/omega_m);  // effectively omega_m -> h .> H0
    H0 = 100*h;

    ochh = omhh_local-obhh_local;
    
    //A_k005  = A_k_WMAP * pow(25.0, n-1); // used by stupid WMAP -> used by CAMB

    /*************************************************/
    /** Write out the initialization parameter file **/
    /*************************************************/
    ifp=fopen("my_params.ini", "w"); 
    write_params_ini(ifp, do_nonlinear, H0, obhh_local, ochh, n, dn_dlnk, A_k005, w0,r_s);
    fclose(ifp); 

    /******************************************/
    /** Run camb; really camb my_params.ini **/
    /*****************************************/
    system("~/new_code/CAMB/camb my_params.ini");
    //system("camb my_params.ini");
    
    /****************************************************************************/
    /** Convert T(k) into properly normalized T(k) which is 1 on large scales  **/
    /**                       Note that k is in h/Mpc                          **/
    /****************************************************************************/

    int N_LINES = file_eof_linecount("transfer_out.dat");

    /****************************************************************************/
    /** Check if the first line starts with hashtag (ie. has variable names)   **/
    /** and skip that first line if that's the case                            **/
    /****************************************************************************/    
    int N_SKIP = nskip_first_line_of_file_or_not("transfer_out.dat");
    
    N_LINES =  N_LINES - N_SKIP;

    /*****************************/
    /** skip the unwanted lines **/
    /*****************************/
    ifp=fopen("transfer_out.dat", "r");
    int N_BUFFER=1e6;
    char unwanted_lines[N_BUFFER]; 
    for(i=0; i < N_SKIP; i++)
        fgets(unwanted_lines, N_BUFFER, ifp);    

    /*********************************************************************************/
    /* if necessary, read extra new columns in damn new CAMB version of T(k) output **/
    /*********************************************************************************/
    int ncols = count_columns_in_last_line_file("transfer_out.dat");

    int j;
    double dum;
    
    ifp1=fopen("my_transfer_out.dat", "w");
    for (i=1; i<=N_LINES; i++)
    {
        fscanf(ifp,"%lf %lf %lf %lf %lf %lf %lf ",
               &k, &Tk_CDM, &Tk_bar, &Tk_phot, &Tk_nu_massless, &Tk_nu_massive, &Tk_tot);
        if (ncols > 7) 
            for (j=8; j<= ncols; j++)
                fscanf(ifp, "%lf ", &dum);
        fscanf(ifp,"\n");
        
        if (i == 1) Tk_norm = Tk_tot;
        fprintf(ifp1, "%e %e\n", k, Tk_tot/Tk_norm);
        //printf("%e %e\n", k, Tk_tot);
    }
    fclose(ifp);
    fclose(ifp1);

    
    /***************/
    /** clean up  **/
    /***************/
    system("rm -rf transfer_out.dat");

}

//Im adding another argument, r_s = redshift to this function. 07/03/2018
void write_params_ini(FILE *ifp, int do_nonlinear, double H0, double obhh_local, double ochh, double ns, 
		      double dn_dlnk, double As, double w, double r_s)
{

      /*Parameters for CAMB */    
      /*output_root is prefixed to output file names */
      fprintf(ifp, "%s\n", "output_root = ");
      
      /*What to do */
      fprintf(ifp, "%s\n", "get_scalar_cls = F");
      fprintf(ifp, "%s\n", "get_vector_cls = F");
      fprintf(ifp, "%s\n", "get_tensor_cls = F");
      fprintf(ifp, "%s\n", "want_CMB = F");

      fprintf(ifp, "%s\n", "get_transfer = T");

      /*#if do_lensing then scalar_output_file contains additional columns of l^4 C_l^{pp} and l^3 C_l^{pT}
        #where p is the projected potential. Output lensed CMB Cls (without tensors) are in lensed_output_file below.*/

      fprintf(ifp, "%s\n", "do_lensing = F");

      /* 0: linear, 1: non-linear matter power (HALOFIT), 2: non-linear CMB lensing (HALOFIT) */
      fprintf(ifp, "%s %d\n", "do_nonlinear =", do_nonlinear);

     /* *#Maximum multipole and k*eta.  */
     /* #  Note that C_ls near l_max are inaccurate (about 5%), go to 50 more than you need */
     /* #  Lensed power spectra are computed to l_max_scalar-250 where accurate at %-level */
     /* #  For high accuracy lensed spectra set l_max_scalar = (l you need) + 500 */
     /* #  To get accurate lensed BB need to have l_max_scalar>2000, k_eta_max_scalar > 10000 */
     /* #  Otherwise k_eta_max_scalar=2*l_max_scalar usually suffices */
      fprintf(ifp, "%s %d\n", "l_max_scalar =", 2000);
      fprintf(ifp, "%s %d\n", "k_eta_max_scalar =", 4000);

      /*#   Tensor settings should be less than or equal to the above */
      fprintf(ifp, "%s %d\n", "l_max_tensor =",  1500);
      fprintf(ifp, "%s %d\n", "k_eta_max_tensor =", 3000);
           
      /*Main cosmological parameters, neutrino masses are assumed degenerate */
      /* If use_phyical set phyiscal densities in baryone, CDM and neutrinos + Omega_k */
      
      fprintf(ifp, "%s\n", "use_physical = T");
      fprintf(ifp, "%s %f\n", "ombh2 =", obhh_local);
      fprintf(ifp, "%s %f\n", "omch2 =", ochh);
      fprintf(ifp, "%s   \n", "omnuh2 = 0.000"); // formerly zero
      fprintf(ifp, "%s   \n", "omk = 0");
      fprintf(ifp, "%s %f\n", "hubble =", H0);

      /* effective equation of state parameter for dark energy, assumed constant */
      fprintf(ifp, "%s %f\n", "w = ", w);
     /* constant comoving sound speed of the dark energy (1=quintessence) */
      fprintf(ifp, "%s\n", "cs2_lam = 1");

      
      /*massless_neutrinos is the effective number (for QED + non-instantaneous decoupling) */
      fprintf(ifp, "%s\n", "temp_cmb = 2.725");
      fprintf(ifp, "%s\n", "helium_fraction = 0.24");
      fprintf(ifp, "%s\n", "massless_neutrinos = 2.04"); // formerly 3.03
      fprintf(ifp, "%s\n", "massive_neutrinos = 1");     // formerly 0
      
      fprintf(ifp, "%s\n", "nu_mass_eigenstates = 1");
      fprintf(ifp, "%s\n", "nu_mass_degeneracies = 1");
      fprintf(ifp, "%s\n", "nu_mass_fractions = 1");
      
      /*Reionization (assumed sharp), ignored unless reionization = T */
      fprintf(ifp, "%s\n", "reionization = F");
      fprintf(ifp, "%s\n", "re_use_optical_depth = F");
      fprintf(ifp, "%s\n", "re_optical_depth     = 0.04");
      /*If re_use_optical_depth = F then use following, otherwise ignored */
      fprintf(ifp, "%s\n", "re_redshift          = 12");
      fprintf(ifp, "%s\n", "re_ionization_frac   = 1");

      
      fprintf(ifp, "%s\n", "pivot_scalar =   5.00000007450580597E-002");
      fprintf(ifp, "%s\n", "pivot_tensor =   5.00000007450580597E-002");
      
      /*Initial power spectrum, amplitude, spectral index and running */
      /* NOTE that CAMB and WMAP amplitudes are same - both defined at k=0.05*/
      fprintf(ifp, "%s\n", "initial_power_num = 1");
      fprintf(ifp, "%s %e\n", "scalar_amp(1) = ", As);
      fprintf(ifp, "%s %f\n", "scalar_spectral_index(1) =", ns);
      fprintf(ifp, "%s %f\n", "scalar_nrun(1) =", dn_dlnk);

      /* ratio is that of the initial tens/scal power spectrum amplitudes*/
      fprintf(ifp, "%s\n", "initial_ratio(1)          = 1");
      /* note vector modes use the scalar settings above*/


      /*Initial scalar perturbation mode (adiabatic=1, CDM iso=2, Baryon iso=3, */
      /* neutrino density iso =4, neutrino velocity iso = 5)  */
      fprintf(ifp, "%s\n", "initial_condition = 1");
      
      /* If above is zero, use modes in the following (totally correlated) proportions */
      /* Note: we assume all modes have the same initial power spectrum*/
      fprintf(ifp, "%s\n", "initial_vector = -1 0 0 0 0");

      /* For vector modes: 0 for regular (neutrino vorticity mode), 1 for magnetic */
      fprintf(ifp, "%s\n", "vector_mode = 0");

      /* Normalization*/
      fprintf(ifp, "%s\n", "COBE_normalize = F");

      /* Transfer function settings, transfer_kmax=0.5 is enough for sigma_8*/
      fprintf(ifp, "%s\n", "transfer_high_precision = F");
      fprintf(ifp, "%s\n", "transfer_kmax           = 100.0");
      fprintf(ifp, "%s\n", "transfer_k_per_logint   = 5");
      fprintf(ifp, "%s\n", "transfer_num_redshifts  = 1");

      //Chao: here a user-specified redshift is added
      fprintf(ifp, "%s %f\n", "transfer_redshift(1)    = ", r_s);
      fprintf(ifp, "%s\n", "transfer_filename(1)    = transfer_out.dat");

      /* Matter power spectrum output against k/h in units of h^3 Mpc^{-3}*/
      fprintf(ifp, "%s\n", "transfer_matterpower(1) = matterpower.dat") ;    
      

      
      fprintf(ifp, "%s\n", "recombination = 1");
      
      
      /*Output files not produced if blank. make camb_fits to use use the FITS setting. */
      fprintf(ifp, "%s\n", "scalar_output_file = scalCls.dat"); 
      fprintf(ifp, "%s\n", "lensed_output_file = lensedCls.dat");
      fprintf(ifp, "%s\n", "tensor_output_file = tensCls.dat");
      fprintf(ifp, "%s\n", "total_output_file  = totCls.dat");

 
      /*if number_of_threads=0 assigned automatically */
      fprintf(ifp, "%s\n", "number_of_threads = 0");
      
      /*Default scalar accuracy is about 0.3% (except lensed BB). */
      /*For 0.1%-level try accuracy_boost=2, l_accuracy_boost=2. */
      
      /*Increase accuracy_boost to decrease time steps, use more k values,  etc. */
      /*Decrease to speed up at cost of worse accuracy. Suggest 0.8 to 3. */
      fprintf(ifp, "%s\n", "accuracy_boost = 1");
      
      /*  Larger to keep more terms in the hierarchy evolution.  */
      fprintf(ifp, "%s\n", "l_accuracy_boost = 1");
      
      /*Increase to use more C_l values for interpolation. */
      /*Increasing a bit will improve the polarization accuracy at l up to 200 - */
      /*interpolation errors may be up to 3% */
      /*Decrease to speed up non-flat models a bit */
      fprintf(ifp, "%s\n", "l_sample_boost = 1");

}
int file_eof_linecount(char *filename)
{
    int NL=0, c;
  FILE *ifp;
  ifp=fopen(filename, "r");

  do      
  {
      c = fgetc(ifp);
      if (c == '\n') NL++;
  }
  while (c != EOF);
  
  // printf("Number of lines = %d\n",  NL);
  fclose(ifp);
  return NL;
}
int count_columns_in_last_line_file(char *filename)
{  // count the number of columns in a file
  // because the first line is sometimes variables, count in the last line just in case
   
    char bash_cmd[256];
    char buffer[1000];
    int ncols;
    FILE *pipe;

    // read the last line in the file since it's definitely just numbers
    sprintf(bash_cmd, "tail -1 %s | wc -w", filename);
    
    // using popen to read char by char
    pipe = popen(bash_cmd, "r");

    if (NULL == pipe) {
        perror("pipe error");
        exit(1);
    } 

    fgets(buffer, sizeof(buffer), pipe);

    pclose(pipe);
    //printf("buffer is %s\n", buffer);

    ncols = atoi (buffer);
    //printf("ncols=%d\n", ncols);
    return ncols;
}
int nskip_first_line_of_file_or_not(char *filename)
{
    //returns 1 if first line starts with a hashtag and should be skipped
    //        returns 0 if not
    char first_entry[256];
    FILE *ifp;
    int nskip;
    
    ifp=fopen(filename, "r");
    fscanf(ifp, "%s ", first_entry); // note no ampersand!
    
    if (first_entry[0] == '#')
        nskip = 1;    
    else nskip = 0;

    fclose(ifp);
    
    return nskip;
}
   
