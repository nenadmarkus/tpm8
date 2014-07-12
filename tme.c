#ifndef ABS
#define ABS(x) ((x)>0?(x):(-(x)))
#endif

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#define SWAP(a, b) (((a) == (b)) || (((a) = (a)^(b)), ((b) = (a)^(b)), ((a) = (a)^(b))))

#define SAVE_TEMPLATE(t, file) fwrite((t), sizeof(int32_t), (t)[0]+1, (file))
#define LOAD_TEMPLATE(t, file) fread(&(t)[0], sizeof(int32_t), 1, (file)), ((t)[0]==0)?0:fread(&(t)[1], sizeof(int32_t), (t)[0], (file))

/*
	
*/

#define _FIXED_POINT_SCALE_ 256
#define _SQR_FIXED_POINT_SCALE_ (_FIXED_POINT_SCALE_*_FIXED_POINT_SCALE_)

void compute_rcs_transformation(int* T, int r, int c, int s)
{
	T[0] = _FIXED_POINT_SCALE_*s; T[1] = 0; T[2] = _SQR_FIXED_POINT_SCALE_*r;
	T[3] = 0; T[4] = _FIXED_POINT_SCALE_*s; T[5] = _SQR_FIXED_POINT_SCALE_*c;
}

/*
	
*/

int bintest(int32_t tcode, int threshold, int* T, uint8_t* pixels, int nrows, int ncols, int ldim)
{
	int r1, c1, r2, c2;

	int8_t* p;

	//
	p = (int8_t*)&tcode;

	//
	r1 = (T[0]*p[0] + T[1]*p[1] + T[2])/_SQR_FIXED_POINT_SCALE_;
	c1 = (T[3]*p[0] + T[4]*p[1] + T[5])/_SQR_FIXED_POINT_SCALE_;

	r2 = (T[0]*p[2] + T[1]*p[3] + T[2])/_SQR_FIXED_POINT_SCALE_;
	c2 = (T[3]*p[2] + T[4]*p[3] + T[5])/_SQR_FIXED_POINT_SCALE_;

	//
#ifndef USE_RGB
	return ( ABS(pixels[r1*ldim+c1]-pixels[r2*ldim+c2]) > threshold );
#else
	return
	(
		(ABS(pixels[r1*ldim+3*c1+0]-pixels[r2*ldim+3*c2+0])>threshold)||
		(ABS(pixels[r1*ldim+3*c1+1]-pixels[r2*ldim+3*c2+1])>threshold)||
		(ABS(pixels[r1*ldim+3*c1+2]-pixels[r2*ldim+3*c2+2])>threshold)
	);
#endif
}

int get_bintest_proximity(int32_t t1, int32_t t2)
{
	int r1, c1, r2, c2;

	int8_t* p1;
	int8_t* p2;

	//
	p1 = (int8_t*)&t1;
	p2 = (int8_t*)&t2;

	//
	r1 = (p1[0]+p1[2])/2;
	c1 = (p1[1]+p1[3])/2;

	r2 = (p2[0]+p2[2])/2;
	c2 = (p2[1]+p2[3])/2;

	//
	return (r1-r2)*(r1-r2) + (c1-c2)*(c1-c2);
}

float get_area_orientation(int r, int c, uint8_t pixels[], int nrows, int ncols, int ldim, int ksize)
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

	for(i=MAX(r-s/2, 0); i<MIN(r+s/2, nrows); ++i)
		for(j=MAX(c-s/2, 0); j<MIN(c+s/2, ncols); ++j)
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
	///T[0] = _FIXED_POINT_SCALE_*s; T[1] = 0; T[2] = _SQR_FIXED_POINT_SCALE_*r;
	///T[3] = 0; T[4] = _FIXED_POINT_SCALE_*s; T[5] = _SQR_FIXED_POINT_SCALE_*c;

	compute_rcs_transformation(T, r, c, s);

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

			int e, elist[1024], ok;

			//
			e = mwcrand()%en;

			//
			o = get_area_orientation(ers[e], ecs[e], pixels, nrows, ncols, ldim, 2);

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
			ok = 1;

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

				///Tp[0] = _FIXED_POINT_SCALE_*s; Tp[1] = 0; Tp[2] = _SQR_FIXED_POINT_SCALE_*(r+mwcrand()%(p/2+1)-(p/2));
				///Tp[3] = 0; Tp[4] = _FIXED_POINT_SCALE_*s; Tp[5] = _SQR_FIXED_POINT_SCALE_*(c+mwcrand()%(p/2+1)-(p/2));

				Tp[0] = T[0]; Tp[1] = T[1]; Tp[2] = _SQR_FIXED_POINT_SCALE_*(r+mwcrand()%(p/2+1)-(p/2));
				Tp[3] = T[3]; Tp[4] = T[4]; Tp[5] = _SQR_FIXED_POINT_SCALE_*(c+mwcrand()%(p/2+1)-(p/2));

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
			int ok;

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
			ok = 1;

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

int match_template_at(int32_t template[], int threshold, int* T, int* pn1, int n0max, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int i, n0;

	int8_t* ptr;

	//
	if(!template[0])
		return 1;

	//
	n0 = 0;

	ptr = (int8_t*)&template[1];

	for(i=0; i<template[0]; ++i)
	{
		if( !bintest(template[i+1], threshold, T, pixels, nrows, ncols, ldim) )
		{
			++n0;

			//
			if(n0 > n0max)
				return 0;
		}
	}

	//
	*pn1 = template[0] - n0;

	//
	return 1;
}

/*
	
*/

typedef struct _tnode
{
	int leaf, tag;

	int32_t* template;

	struct _tnode* subtree1;
	struct _tnode* subtree2;

} tnode;

int save_tree_to_file(tnode* root, FILE* file)
{
	int32_t dummy;

	fwrite(&root->leaf, sizeof(int), 1, file);
	fwrite(&root->tag, sizeof(int), 1, file);

	if(root->template)
		SAVE_TEMPLATE(root->template, file);
	else
	{
		dummy = 0;
		fwrite(&dummy, sizeof(int32_t), 1, file);
	}

	if(!root->leaf)
	{
		return 
			save_tree_to_file(root->subtree1, file)
				|
			save_tree_to_file(root->subtree2, file);
	}
	else
		return 1;
}

tnode* load_tree_from_file(FILE* file)
{
	tnode* root = (tnode*)malloc(sizeof(tnode));

	//
	fread(&root->leaf, sizeof(int), 1, file);
	fread(&root->tag, sizeof(int), 1, file);

	//
	root->template = (int32_t*)malloc((MAXNUMTESTS+1)*sizeof(int32_t));

	LOAD_TEMPLATE(root->template, file);

	//
	if(!root->leaf)
	{
		root->subtree1 = load_tree_from_file(file);
		root->subtree2 = load_tree_from_file(file);
	}

	//
	return root;
}

int learn_joint_features(int32_t stack[], int stacksize, int maxstacksize, int* Ts[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n, int32_t tcodepool[], int tcodepoolsize, int perturbationstrength)
{
	int i, j, k, newstacksize, numiters, maxnumiters;

	//
	newstacksize = stacksize;

	maxnumiters = 2048;

	numiters = 0;

	while(newstacksize<maxstacksize && numiters<maxnumiters)
	{
		int ok, nfails;

		//
		k = mwcrand()%tcodepoolsize;

		stack[newstacksize] = tcodepool[k];

		//
		ok = 1;

		for(i=0; i<newstacksize; ++i)
			/*
				proximity requirements
			*/
			if(get_bintest_proximity(stack[i], stack[newstacksize]) < (maxnumiters-numiters)*255/maxnumiters)
				ok = 0;

		nfails = 0;

		for(k=0; k<n; ++k)
		{
			/*
				stability requirements
			*/

			int* T = Ts[k];

			//
			if(0==bintest(stack[newstacksize], THRESHOLD, T, pixelss[inds[k]], nrowss[inds[k]], ncolss[inds[k]], ldims[inds[k]]))
				nfails += 1000;

			for(i=0; i<32; ++i)
			{
				int Tp[6];

				Tp[0] = T[0]; Tp[1] = T[1]; Tp[2] = T[2] + _SQR_FIXED_POINT_SCALE_*( mwcrand()%(perturbationstrength+1) - perturbationstrength );
				Tp[3] = T[3]; Tp[4] = T[4]; Tp[5] = T[5] + _SQR_FIXED_POINT_SCALE_*( mwcrand()%(perturbationstrength+1) - perturbationstrength );

				if( 0==bintest(stack[newstacksize], THRESHOLD, Tp, pixelss[inds[k]], nrowss[inds[k]], ncolss[inds[k]], ldims[inds[k]]) )
				{
					++nfails;
					break;
				}
			}
		}

		if(nfails > 0)
			ok = 0;

		//
		if(ok)
			++newstacksize;

		//
		++numiters;
	}

	//
	return newstacksize;
}

float get_similarity
	(
		int32_t t1[], int* T1, uint8_t p1[], int nrows1, int ncols1, int ldim1,
		int32_t t2[], int* T2, uint8_t p2[], int nrows2, int ncols2, int ldim2
	)
{
	int n1, n2;
	int s12, s21;

	//
	n1 = t1[0];
	n2 = t2[0];

	//
	match_template_at(t1, THRESHOLD, T2, &s12, n1, p2, nrows2, ncols2, ldim2);
	match_template_at(t2, THRESHOLD, T1, &s21, n2, p1, nrows1, ncols1, ldim1);

	//
	return ( s12/(float)n1 + s21/(float)n2 )/2.0f;
}

int partition_data(int32_t* templates[], int* Ts[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n)
{
	int stop;
	int i, j;

	//
	if(n == 2)
		return 1;

	//
	stop = 0;

	i = 0;
	j = n - 1;

	while(!stop)
	{
		//
		while
		(
			get_similarity
			(
				templates[inds[0]], Ts[inds[0]], pixelss[inds[0]], nrowss[inds[0]], ncolss[inds[0]], ldims[inds[0]],
				templates[inds[i]], Ts[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]
			)
				>=
			get_similarity
			(
				templates[inds[i]], Ts[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]],
				templates[inds[n-1]], Ts[inds[n-1]], pixelss[inds[n-1]], nrowss[inds[n-1]], ncolss[inds[n-1]], ldims[inds[n-1]]
			)
		)
		{
			if( i==j )
				break;
			else
				++i;
		}

		while
		(
			get_similarity
			(
				templates[inds[0]], Ts[inds[0]], pixelss[inds[0]], nrowss[inds[0]], ncolss[inds[0]], ldims[inds[0]],
				templates[inds[j]], Ts[inds[j]], pixelss[inds[j]], nrowss[inds[j]], ncolss[inds[j]], ldims[inds[j]]
			)
				<=
			get_similarity
			(
				templates[inds[j]], Ts[inds[j]], pixelss[inds[j]], nrowss[inds[j]], ncolss[inds[j]], ldims[inds[j]],
				templates[inds[n-1]], Ts[inds[n-1]], pixelss[inds[n-1]], nrowss[inds[n-1]], ncolss[inds[n-1]], ldims[inds[n-1]]
			)
		)
		{
			if( i==j )
				break;
			else
				--j;
		}

		//
		if( i==j )
			stop = 1;
		else
			SWAP(inds[i], inds[j]);
	}

	//
	return i; // ?
}

tnode* grow_subtree(int depth, int32_t stack[], int stacksize, int maxnumtests, int* Ts[], int32_t* templates[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n, int32_t tcodepool[], int tcodepoolsize, int perturbationstrength)
{
	int i, newstacksize, n1, n2;

	tnode* root = 0;

	//
	if(n==0)
		return 0;

	root = (tnode*)malloc(sizeof(tnode));

	//
	newstacksize = learn_joint_features(stack, stacksize, stacksize+maxnumtests, Ts, pixelss, nrowss, ncolss, ldims, inds, n, tcodepool, tcodepoolsize, perturbationstrength);

	if(newstacksize-stacksize > maxnumtests/2)
	{
		root->template = (int32_t*)malloc((maxnumtests+1)*sizeof(int32_t));

		//
		root->template[0] = newstacksize - stacksize;

		for(i=stacksize; i<newstacksize; ++i)
			root->template[1+i-stacksize] = stack[i];

		//
		///if(n>1)printf("%d: %d, %d\n", depth, root->template[0], n);
	}
	else
		root->template = 0;

	if(n == 1)
	{
		root->leaf = 1;
		root->tag = inds[0];

		root->subtree1 = 0;
		root->subtree2 = 0;

		//printf("%d %d\n", depth, newstacksize);

		return root;
	}
	else
	{
		root->leaf = 0;
		root->tag = 0;
	}

	// split data
	i = mwcrand()%n;
	SWAP(inds[0], inds[i]);

	i = 1+mwcrand()%(n-1);
	SWAP(inds[n-1], inds[i]);

	n1 = partition_data(templates, Ts, pixelss, nrowss, ncolss, ldims, inds, n); // 0th and (n-1)st samples serve as "anchors"

	n1 = MAX(n1, 1); // hack?

	n2 = n - n1;

	//
	root->subtree1 = grow_subtree(depth+1, stack, newstacksize, maxnumtests, Ts, templates, pixelss, nrowss, ncolss, ldims, &inds[0 ], n1, tcodepool, tcodepoolsize, perturbationstrength);
	root->subtree2 = grow_subtree(depth+1, stack, newstacksize, maxnumtests, Ts, templates, pixelss, nrowss, ncolss, ldims, &inds[n1], n2, tcodepool, tcodepoolsize, perturbationstrength);

	//
	return root;
}

tnode* grow_tree(int* Ts[], int32_t* templates[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int n, int32_t tcodepool[], int tcodepoolsize, int perturbationstrength)
{
	int i;
	int* inds;

	int maxstacksize;
	int32_t* stack;

	tnode* root;

	//
	inds = (int*)malloc(n*sizeof(int));

	for(i=0; i<n; ++i)
		inds[i] = i;

	//
	maxstacksize = n*MAXNUMTESTS;

	stack = (int32_t*)malloc(maxstacksize*sizeof(int32_t));

	//
	root = grow_subtree(0, stack, 0, MAXNUMTESTS/4, Ts, templates, pixelss, nrowss, ncolss, ldims, inds, n, tcodepool, tcodepoolsize, perturbationstrength);

	//
	free(stack);

	//
	return root;
}

int get_tree_output(tnode* root, int threshold, int n0max, int* T, uint8_t pixels[], int nrows, int ncols, int ldim, int tags[], int* n)
{
	if(root->template)
	{
		int n1;

		if(!match_template_at(root->template, threshold, T, &n1, n0max, pixels, nrows, ncols, ldim))
			return 0;
	}

	if(root->leaf)
	{
		tags[*n] = root->tag;
		++*n;

		return 1;
	}
	else
	{
		return
			get_tree_output(root->subtree1, threshold, n0max, T, pixels, nrows, ncols, ldim, tags, n)
				|
			get_tree_output(root->subtree2, threshold, n0max, T, pixels, nrows, ncols, ldim, tags, n);
	}
}