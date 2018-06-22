void my_locate(double xx[], unsigned long n, double x, unsigned long *j)
{
	unsigned long ju,jm,jl;
	int ascnd;

	jl=0;
	ju=n+1;
	ascnd=(xx[n] >= xx[1]);
	while (ju-jl > 1) {
		jm=(ju+jl) >> 1;
		if (x >= xx[jm] == ascnd)
			jl=jm;
		else
			ju=jm;
	}
	if (x == xx[1]) *j=1;
	else if(x == xx[n]) *j=n-1;
	else *j=jl;
}
