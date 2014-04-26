#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#define SAVE_TEMPLATE(t, file) fwrite((t), sizeof(int32_t), (t)[0]+1, (file))
#define LOAD_TEMPLATE(t, file) fread(&(t)[0], sizeof(int32_t), 1, (file)), ((t)[0]==0)?0:fread(&(t)[1], sizeof(int32_t), (t)[0], (file))

float getorient(int r, int c, uint8_t pixels[], int nrows, int ncols, int ldim, int ksize)
{
	int r1, r2, c1, c2, cnt;
	float norm, vdiff;

	float gr, gc, o;

	//
	cnt = 0;

	gr = 0.0f;
	gc = 0.0f;

	for(r1=-ksize; r1<+ksize; ++r1)
		for(c1=-ksize; c1<+ksize; ++c1)
		{
			if(r1*r1+c1*c1 > ksize*ksize)
				continue;

			for(r2=-ksize; r2<+ksize; ++r2)
				for(c2=-ksize; c2<+ksize; ++c2)
				{
					if(r2*r2+c2*c2 > ksize*ksize)
						continue;

					norm = sqrt( (r1-r2)*(r1-r2) + (c1-c2)*(c1-c2) );

					if(norm < 3)
						continue;

					vdiff = pixels[(r+r1)*ldim+(c+c1)] - pixels[(r+r2)*ldim+(c+c2)];
				
					gr += vdiff*(r1-r2)/norm;
					gc += vdiff*(c1-c2)/norm;

					++cnt;
				}
		}

	//
	gr /= cnt;
	gc /= cnt;

	//
	if(gc == 0.0f)
		return 1.57f;
	else
		return atan( gr/gc );
}

int learn_template(int32_t template[], int maxnumtests, int useorientation, float s2p, int r, int c, int s, uint8_t pixels[], uint8_t mask[], int nrows, int ncols, int ldim, int threshold)
{
	int i, j, n, numiters, maxnumiters, p;

	int8_t* ptr;
	int32_t b;

	static int ers[400*400], ecs[400*400], en;

	int r1, c1, r2, c2;

	int T[6];

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

	// compute the transformation matrix
	T[0] = _FIXED_POINT_SCALE_*s; T[1] = 0; T[2] = _SQR_FIXED_POINT_SCALE_*r;
	T[3] = 0; T[4] = _FIXED_POINT_SCALE_*s; T[5] = _SQR_FIXED_POINT_SCALE_*c;

	//
	n = 0;
	ptr = (int8_t*)&template[1];

	if(useorientation == 1)
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
			o = getorient(ers[e], ecs[e], pixels, nrows, ncols, ldim, 2);

			//
			r1 = MIN(MAX(r-s/2+1, ers[e]-sin(o)*p), r+s/2-1);
			c1 = MIN(MAX(c-s/2+1, ecs[e]-cos(o)*p), c+s/2-1);

			r2 = MIN(MAX(r-s/2+1, ers[e]+sin(o)*p), r+s/2-1);
			c2 = MIN(MAX(c-s/2+1, ecs[e]+cos(o)*p), c+s/2-1);

			//
			ptr[4*n+0] = _FIXED_POINT_SCALE_*(r1-r)/s;
			ptr[4*n+1] = _FIXED_POINT_SCALE_*(c1-c)/s;
			ptr[4*n+2] = _FIXED_POINT_SCALE_*(r2-r)/s;
			ptr[4*n+3] = _FIXED_POINT_SCALE_*(c2-c)/s;

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

			if( 0==bintest(b, threshold, T, pixels, nrows, ncols, ldim) )
				ok = 0;

			for(i=0; i<32; ++i)
			{
				/*
					stability requirements
				*/
				int Tp[6];

				Tp[0] = _FIXED_POINT_SCALE_*s; Tp[1] = 0; Tp[2] = _SQR_FIXED_POINT_SCALE_*(r+mwcrand()%(p/2+1)-(p/2));
				Tp[3] = 0; Tp[4] = _FIXED_POINT_SCALE_*s; Tp[5] = _SQR_FIXED_POINT_SCALE_*(c+mwcrand()%(p/2+1)-(p/2));

				if( 0==bintest(b, threshold, Tp, pixels, nrows, ncols, ldim) )
					ok = 0;
			}

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
			ptr[4*n+0] = _FIXED_POINT_SCALE_*(r1-r)/s;
			ptr[4*n+1] = _FIXED_POINT_SCALE_*(c1-c)/s;
			ptr[4*n+2] = _FIXED_POINT_SCALE_*(r2-r)/s;
			ptr[4*n+3] = _FIXED_POINT_SCALE_*(c2-c)/s;

			b = *(int32_t*)&ptr[4*n+0];

			//
			int ok = 1;

			for(i=0; i<32; ++i)
			{
				int Tp[6];

				Tp[0] = _FIXED_POINT_SCALE_*s; Tp[1] = 0; Tp[2] = _SQR_FIXED_POINT_SCALE_*(r-p+mwcrand()%(2*p));
				Tp[3] = 0; Tp[4] = _FIXED_POINT_SCALE_*s; Tp[5] = _SQR_FIXED_POINT_SCALE_*(c-p+mwcrand()%(2*p));

				if( 1==bintest(b, threshold, Tp, pixels, nrows, ncols, ldim) )
					ok = 0;
			}

			if(ok)
				++n;

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

	int T[6], norm;

	//
	if(!template[0])
		return 1;

	// compute the transformation matrix
	T[0] = _FIXED_POINT_SCALE_*s; T[1] = 0; T[2] = _SQR_FIXED_POINT_SCALE_*r;
	T[3] = 0; T[4] = _FIXED_POINT_SCALE_*s; T[5] = _SQR_FIXED_POINT_SCALE_*c;

	//
	n0 = 0;
	r0 = 0;

	ptr = (int8_t*)&template[1];

	for(i=0; i<template[0]; ++i)
	{
		if( !bintest(template[i+1], threshold, T, pixels, nrows, ncols, ldim) )
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