CC = gcc
CFLAGS=-g -Wall -fopenmp
IDIR   = -I/usr/common/software/gsl/2.1/intel/include      
LDIR   = -L/usr/common/software/gsl/2.1/intel/lib      -L ./recipes/lib
cl:  cl_theory.o functions.o run_camb_get_Tk.o  spline_functions.o get_sn_cov.o
	$(CC) $(CFLAGS) -o  cl cl_theory.o functions.o run_camb_get_Tk.o spline_functions.o get_sn_cov.o $(LDIR) -lgsl -lgslcblas -lrecipes -lm
sig:   signal.o functions.o run_camb_get_Tk.o  spline_functions.o get_sn_cov.o 
	$(CC) $(CFLAGS) -o  sig signal.o functions.o run_camb_get_Tk.o spline_functions.o get_sn_cov.o $(LDIR) -lgsl -lgslcblas -lrecipes -lm
galcov:  galcov_main.o functions.o run_camb_get_Tk.o  spline_functions.o get_gal_cov.o get_sn_cov.o
	$(CC) $(CFLAGS) -o  galcov galcov_main.o functions.o run_camb_get_Tk.o spline_functions.o get_sn_cov.o get_gal_cov.o $(LDIR) -lgsl -lgslcblas -lrecipes -lm
fish:   fisher_est_fs8.o functions.o run_camb_get_Tk.o spline_functions.o get_sn_cov.o 
	$(CC) $(CFLAGS) -o  fish fisher_est_fs8.o functions.o run_camb_get_Tk.o  spline_functions.o  get_sn_cov.o $(LDIR) -lgsl -lgslcblas -lrecipes  -lchealpix -lm 
fish:   fisher_est_fs8.o functions.o run_camb_get_Tk.o spline_functions.o get_sn_cov.o 
	$(CC) $(CFLAGS) -o  fish fisher_est_fs8.o functions.o run_camb_get_Tk.o  spline_functions.o  get_sn_cov.o $(LDIR) -lgsl -lgslcblas -lrecipes  -lchealpix -lm
test_NR_2d_spline  :  test_NR_2d_spline.o
	$(CC) $(CFLAGS) -o  test_NR_2d_spline test_NR_2d_spline.o $(LDIR) -lrecipes -lm
find_2D_contour: find_2D_contour.o
	$(CC) $(CFLAGS) -o  find_2D_contour find_2D_contour.o $(LDIR) -lrecipes -lm
.c.o:
	$(CC) -c -o $@ $(CFLAGS) $(IDIR) $<
clean:
	-rm -f *.o
