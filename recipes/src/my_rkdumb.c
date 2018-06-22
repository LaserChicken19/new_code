#include "nrutil.h"

double **my_y,*my_xx;

void my_rkdumb(double vstart[], int nvar, double x1, double x2, int nstep,
	void (*my_derivs)(double, double [], double []))
{
        void my_rk4(double my_y[], double dydx[], int n, double x, double h, 
		 double yout[], void (*my_derivs)(double, double [], double []));
	int i,k;
	double x,h;
	double *v,*vout,*dv;

	v=dvector(1,nvar);
	vout=dvector(1,nvar);
	dv=dvector(1,nvar);
	for (i=1;i<=nvar;i++) {
		v[i]=vstart[i];
		my_y[i][1]=v[i];
	}
	my_xx[1]=x1;
	x=x1;
	h=(x2-x1)/nstep;
	for (k=1;k<=nstep;k++) {
		(*my_derivs)(x,v,dv);
		my_rk4(v,dv,nvar,x,h,vout,my_derivs);
		if ((double)(x+h) == x) nrerror("Step size too small in routine rkdumb");
		x += h;
		my_xx[k+1]=x;
		for (i=1;i<=nvar;i++) {
			v[i]=vout[i];
			my_y[i][k+1]=v[i];
		}
	}
	free_dvector(dv,1,nvar);
	free_dvector(vout,1,nvar);
	free_dvector(v,1,nvar);
}
