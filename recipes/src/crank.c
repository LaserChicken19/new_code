void crank(unsigned long n, float w[], float *s)
{
/*
 *  The "unsigned" declaration causes a problem with cc on vax
 */
#ifdef vax
        long j=1,ji,jt;
#else
	unsigned long j=1,ji,jt;
#endif
	float t,rank;

	*s=0.0;
	while (j < n) {
		if (w[j+1] != w[j]) {
			w[j]=j;
			++j;
		} else {
			for (jt=j+1;jt<=n && w[jt]==w[j];jt++);
			rank=0.5*(j+jt-1);
			for (ji=j;ji<=(jt-1);ji++) w[ji]=rank;
			t=jt-j;
			*s += t*t*t-t;
			j=jt;
		}
	}
	if (j == n) w[n]=n;
}
