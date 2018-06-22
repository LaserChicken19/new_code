#include <math.h>
#include "nrutil.h"
#define ITMAX 200

void my_powell(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	double (*func)(double []))
{
        void my_linmin(double p[], double xi[], int n, double *fret,
		    double (*func)(double []));
	int i,ibig,j;
	double del,fp,fptt,t,*pt,*ptt,*xit;

	pt=dvector(1,n);
	ptt=dvector(1,n);
	xit=dvector(1,n);
	*fret=(*func)(p);
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			my_linmin(p,xit,n,fret,func);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			free_dvector(xit,1,n);
			free_dvector(ptt,1,n);
			free_dvector(pt,1,n);
			return;
		}
		if (*iter == ITMAX) 
		  {
		    printf("my_powell exceeds maximum iterations.\n");
		    /* set the function "minimum" value to something huge */
		    *fret=1.0e50;
		    free_dvector(xit,1,n);
		    free_dvector(ptt,1,n);
		    free_dvector(pt,1,n);
		    return;
		  }
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				my_linmin(p,xit,n,fret,func);
				for (j=1;j<=n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}
#undef ITMAX
