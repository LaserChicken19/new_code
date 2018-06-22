#include <math.h>
#include <stdio.h>

void my_choldc(double **a, int n, double p[])
{
	void nrerror(char error_text[]);
	int i,j,k;
	double sum;

	for (i=1;i<=n;i++) {
		for (j=i;j<=n;j++) {
			for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
			if (i == j) {
				if (sum <= 0.0)
				  {
					printf("in choldc sum=%e\n", sum);
					printf("at i=j=%d\n", i);
					printf("a[i][j]=%e\n", a[i][j]);
			for (sum=a[i][j],k=i-1;k>=1;k--) 
				printf("k=%d a[i][k]=%e a[j][k]=%e\n", k,
					a[i][k], a[j][k]);
					nrerror("choldc failed");
				  }
				p[i]=sqrt(sum);
			} else a[j][i]=sum/p[i];
		}
	}
}
