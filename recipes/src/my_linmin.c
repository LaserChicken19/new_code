#include "nrutil.h"
#define TOL 2.0e-4

int ncom;
double *pcom,*xicom,(*nrfunc)(double []);

void my_linmin(double p[], double xi[], int n, double *fret,
	double (*func)(double []))
{
        double my_brent(double ax, double bx, double cx,
		    double (*f)(double), double tol, double *xmin);
	double my_f1dim(double x);
	void my_mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb,
		    double *fc, double (*func)(double));
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=dvector(1,n);
	xicom=dvector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	my_mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,my_f1dim);
	*fret=my_brent(ax,xx,bx,my_f1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_dvector(xicom,1,n);
	free_dvector(pcom,1,n);
}
#undef TOL
