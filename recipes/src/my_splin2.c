#include "nrutil.h"

void my_splin2(double x1a[], double x2a[], double **ya, double **y2a, int m, int n,
	       double x1, double x2, double *y)
{
        void my_spline(double x[], double y[], int n, double yp1, 
		       double ypn, double y2[]);
	void my_splint(double xa[], double ya[], double y2a[], int n, 
		       double  x, double *y);
	int j;
	double *ytmp,*yytmp;

	ytmp=dvector(1,m);
	yytmp=dvector(1,m);
	for (j=1;j<=m;j++)
		my_splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
	my_spline(x1a,yytmp,m,1.0e30,1.0e30,ytmp);
	my_splint(x1a,yytmp,ytmp,m,x1,y);
	free_dvector(yytmp,1,m);
	free_dvector(ytmp,1,m);
}
