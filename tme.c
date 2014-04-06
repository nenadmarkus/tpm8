#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#ifndef ABS
#define ABS(x) ((x)>0?(x):(-(x)))
#endif

#define NORMALIZATION ( 255 )

int bintest(int32_t tcode, int threshold, int r, int c, int s, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int r1, c1, r2, c2;

	int8_t* ptr;

	//
	ptr = (int8_t*)&tcode;

	//
	r1 = (NORMALIZATION*r + ptr[0]*s)/NORMALIZATION;
	c1 = (NORMALIZATION*c + ptr[1]*s)/NORMALIZATION;

	r2 = (NORMALIZATION*r + ptr[2]*s)/NORMALIZATION;
	c2 = (NORMALIZATION*c + ptr[3]*s)/NORMALIZATION;

	///printf("%d %d %d %d\n", r1, c1, r2, c2); getchar();

	//
	return ( ABS(pixels[r1*ldim+c1]-pixels[r2*ldim+c2]) > threshold );
}

int learn_template(int32_t template[], int maxntests, float s2p, int r, int c, int s, uint8_t pixels[], uint8_t edgemap[], int nrows, int ncols, int ldim, int threshold)
{
	int i, j;

	int8_t* ptr;

	int ers[8192], ecs[8192], en;

	//
	en = 0;

	for(i=r-s/2; i<r+s/2; ++i)
		for(j=c-s/2; j<c+s/2; ++j)
			if(edgemap[i*ldim+j])
			{
				if( (i-r)*(i-r) + (j-c)*(j-c) < (s-2)*(s-2)/4 )
				{
					ers[en] = i;
					ecs[en] = j;

					++en;
				}
			}

	//
	int n, niters, maxniters;

	maxniters = 16*maxntests;

	n = 0;
	niters = 0;

	ptr = (int8_t*)&template[1];

	if(en)
	{
		int p = (int)( s*s2p );

		while(n<maxntests && niters<maxniters)
		{
			float o, gr, gc;

			int e;
			int elist[1024];

			//
			e = mwcrand()%en;

			//
			gr = pixels[(ers[e]+1)*ldim+ecs[e]] - pixels[(ers[e]-1)*ldim+ecs[e]];
			gc = pixels[ers[e]*ldim+(ecs[e]+1)] - pixels[ers[e]*ldim+(ecs[e]-1)];

			if(gc == 0)
				o = 1.57f;
			else
				o = atan( gr/gc );

			//
			int r1, c1, r2, c2;

			r1 = ers[e] - sin(o)*p;
			c1 = ecs[e] - cos(o)*p;

			r2 = ers[e] + sin(o)*p;
			c2 = ecs[e] + cos(o)*p;

			//
			r1 = MIN(MAX(r-s/2+1, r1), r+s/2-1);
			c1 = MIN(MAX(c-s/2+1, c1), c+s/2-1);

			r2 = MIN(MAX(r-s/2+1, r2), r+s/2-1);
			c2 = MIN(MAX(c-s/2+1, c2), c+s/2-1);

			//
			if(ABS(pixels[r1*ldim+c1]-pixels[r2*ldim+c2]) > threshold)
			{
				int ok;

				//
				ok = 1;

				// proximity requirements
				for(i=0; i<n; ++i)
				{
					int d2 = (ers[elist[i]]-ers[e])*(ers[elist[i]]-ers[e]) + (ecs[elist[i]]-ecs[e])*(ecs[elist[i]]-ecs[e]);

					if(d2 < (maxniters-niters)*s/maxniters)
						ok = 0;
				}

				// stability requirements
				int nok = 0;

				for(i=0; i<32; ++i)
				{
					int pr1, pc1, pr2, pc2;

					//
					pr1 = r1 + mwcrand()%(p/3+1)-(p/3);
					pc1 = c1 + mwcrand()%(p/3+1)-(p/3);

					pr2 = r2 + mwcrand()%(p/3+1)-(p/3);
					pc2 = c2 + mwcrand()%(p/3+1)-(p/3);

					//
					if(ABS(pixels[pr1*ldim+pc1]-pixels[pr2*ldim+pc2]) > threshold)
						++nok;
				}

				if(nok/32.0f < 0.95f)
					ok = 0;

				//
				if(ok)
				{
					//
					ptr[4*n+0] = NORMALIZATION*(r1-r)/s;
					ptr[4*n+1] = NORMALIZATION*(c1-c)/s;
					ptr[4*n+2] = NORMALIZATION*(r2-r)/s;
					ptr[4*n+3] = NORMALIZATION*(c2-c)/s;

					elist[n] = e;

					//
					++n;
				}
			}

			//
			++niters;
		}
	}

	template[0] = n;

	//
	return n;
}

int match_template_at(int32_t template[], int threshold, int r, int c, int s, int* pn1, int n0max, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int i, n0, r0;

	int8_t* ptr;

	//
	if(!template[0])
		return 0;

	//
	n0 = 0;
	r0 = 0;

	ptr = (int8_t*)&template[1];

	for(i=0; i<template[0]; ++i)
	{
		if( !bintest(template[i+1], threshold, r, c, s, pixels, nrows, ncols, ldim) )
		{
			++n0;
			++r0;

			//
			if(n0 > n0max)
				return 0;

			//
			if(r0 > 3)
				return 0;
		}
		else
			r0 = MAX(r0-1, 0);
	}

	//
	*pn1 = template[0] - n0;

	//
	return 1;
}