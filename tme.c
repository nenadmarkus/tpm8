#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#define SAVE_TEMPLATE(t, file) fwrite((t), sizeof(int32_t), (t)[0]+1, (file));
#define LOAD_TEMPLATE(t, file) fread(&(t)[0], sizeof(int32_t), 1, (file)), fread(&(t)[1], sizeof(int32_t), (t)[0], (file));

int learn_template(int32_t template[], int maxnumtests, float s2p, int r, int c, int s, uint8_t pixels[], uint8_t edgemap[], int nrows, int ncols, int ldim, int threshold)
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
	int n, numiters, maxnumiters;

	maxnumiters = 16*maxnumtests;

	n = 0;
	numiters = 0;

	ptr = (int8_t*)&template[1];

	if(en)
	{
		int p = (int)( s*s2p );

		while(n<maxnumtests && numiters<maxnumiters)
		{
			float o, gr, gc;

			int e;
			int elist[1024];

			int32_t b;

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

			r1 = MIN(MAX(r-s/2+1, ers[e]-sin(o)*p), r+s/2-1);
			c1 = MIN(MAX(c-s/2+1, ecs[e]-cos(o)*p), c+s/2-1);

			r2 = MIN(MAX(r-s/2+1, ers[e]+sin(o)*p), r+s/2-1);
			c2 = MIN(MAX(c-s/2+1, ecs[e]+cos(o)*p), c+s/2-1);

			//
			ptr[4*n+0] = NORMALIZATION*(r1-r)/s;
			ptr[4*n+1] = NORMALIZATION*(c1-c)/s;
			ptr[4*n+2] = NORMALIZATION*(r2-r)/s;
			ptr[4*n+3] = NORMALIZATION*(c2-c)/s;

			b = *(int32_t*)&ptr[4*n+0];
			b = b | 0x1;
			*(int32_t*)&ptr[4*n+0] = b;

			//
			int ok;

			//
			ok = 1;

			// proximity requirements
			for(i=0; i<n; ++i)
			{
				int d2 = (ers[elist[i]]-ers[e])*(ers[elist[i]]-ers[e]) + (ecs[elist[i]]-ecs[e])*(ecs[elist[i]]-ecs[e]);

				if(d2 < (maxnumiters-numiters)*s/maxnumiters)
					ok = 0;
			}

			// stability requirements
			for(i=0; i<32; ++i)
				if( !bintest(b, threshold, r+mwcrand()%(p/2+1)-(p/2), c+mwcrand()%(p/2+1)-(p/2), s, pixels, nrows, ncols, ldim) )
					ok = 0;

			//
			if(ok)
			{
				//
				elist[n] = e;

				//
				++n;
			}

			//
			++numiters;
		}
	}

	template[0] = n;

	//
	return n;
}

int match_template_at(int32_t template[], int threshold, int r, int c, int s, int* pn1, int n0max, int r0max, uint8_t pixels[], int nrows, int ncols, int ldim)
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
			if(r0 > r0max)
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