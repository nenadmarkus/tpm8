#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#define SAVE_TEMPLATE(t, file) fwrite((t), sizeof(int32_t), (t)[0]+1, (file));
#define LOAD_TEMPLATE(t, file) fread(&(t)[0], sizeof(int32_t), 1, (file)), fread(&(t)[1], sizeof(int32_t), (t)[0], (file));

int learn_template(int32_t template[], int maxnumtests, int type, float s2p, int r, int c, int s, uint8_t pixels[], uint8_t mask[], int nrows, int ncols, int ldim, int threshold)
{
	int i, j, n, numiters, maxnumiters, p;

	int8_t* ptr;
	int32_t b;

	static int ers[400*400], ecs[400*400], en;

	int r1, c1, r2, c2;

	//
	en = 0;

	for(i=r-s/2; i<r+s/2; ++i)
		for(j=c-s/2; j<c+s/2; ++j)
			if(mask[i*ldim+j])
			{
				if( (i-r)*(i-r) + (j-c)*(j-c) < (s-2)*(s-2)/4 )
				{
					ers[en] = i;
					ecs[en] = j;

					++en;
				}
			}

	if(!en)
	{
		template[0] = 0;
		return 0;
	}

	//
	n = 0;
	ptr = (int8_t*)&template[1];

	if(type == 1)
	{
		p = (int)( s*s2p );
		maxnumiters = 16*maxnumtests;

		numiters = 0;

		while(n<maxnumtests && numiters<maxnumiters)
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

			//
			int ok = 1;

			for(i=0; i<n; ++i)
			{
				/*
					proximity requirements
				*/
				int d2 = (ers[elist[i]]-ers[e])*(ers[elist[i]]-ers[e]) + (ecs[elist[i]]-ecs[e])*(ecs[elist[i]]-ecs[e]);
				if(d2 < (maxnumiters-numiters)*s/maxnumiters)
					ok = 0;
			}

			for(i=0; i<32; ++i)
				/*
					stability requirements
				*/
				if( 0==bintest(b, threshold, r+mwcrand()%(p/2+1)-(p/2), c+mwcrand()%(p/2+1)-(p/2), s, pixels, nrows, ncols, ldim) )
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
	else
	{
		p = (int)( s*s2p );
		maxnumiters = 64*maxnumtests;

		numiters = 0;

		while(n<maxnumtests && numiters<maxnumiters)
		{
			//
			r1 = r - s/2 + mwcrand()%s;
			c1 = c - s/2 + mwcrand()%s;

			r2 = r - s/2 + mwcrand()%s;
			c2 = c - s/2 + mwcrand()%s;

			if(mask[r1*ldim+c1]==0 || mask[r2*ldim+c2]==0)
				continue;

			if( (r1-r2)*(r1-r2)+(c1-c2)*(c1-c2)>s*s/100 || (r1-r2)*(r1-r2)+(c1-c2)*(c1-c2)<s*s/225 )
				continue;

			//
			ptr[4*n+0] = NORMALIZATION*(r1-r)/s;
			ptr[4*n+1] = NORMALIZATION*(c1-c)/s;
			ptr[4*n+2] = NORMALIZATION*(r2-r)/s;
			ptr[4*n+3] = NORMALIZATION*(c2-c)/s;

			b = *(int32_t*)&ptr[4*n+0];

			//
			int ok = 1;

			for(i=0; i<32; ++i)
				if( 1==bintest(b, threshold, r-p+mwcrand()%(2*p), c-p+mwcrand()%(2*p), s, pixels, nrows, ncols, ldim) )
					ok = 0;

			if(ok)
			{
				++n;
			}

			//
			++numiters;
		}
	}

	//
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