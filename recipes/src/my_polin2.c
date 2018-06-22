#include "nrutil.h"

void my_polin2(double x1a[], double x2a[], double **ya, int m, int n,
	      double x1, double x2, double *y, double *dy)
{
        void my_polint(double xa[], double ya[], int n, double x, 
		       double *y, double *dy);
	int j;
	double *ymtmp;

	ymtmp=dvector(1,m);
	for (j=1;j<=m;j++) {
		my_polint(x2a,ya[j],n,x2,&ymtmp[j],dy);
	}
	my_polint(x1a,ymtmp,m,x1,y,dy);
	free_dvector(ymtmp,1,m);
}
