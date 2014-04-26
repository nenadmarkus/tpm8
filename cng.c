#define SWAP(a, b) (((a) == (b)) || (((a) = (a)^(b)), ((b) = (a)^(b)), ((a) = (a)^(b))))

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

int learn_cluster_features(int32_t stack[], int stacksize, int maxstacksize, float s2p, int rs[], int cs[], int ss[], uint8_t* pixelss[], uint8_t* edgemaps[], int nrowss[], int ncolss[], int ldims[], int inds[], int n, int threshold)
{
	int i, j, k, newstacksize, numiters, maxnumiters, p;

	int8_t* stackbyteptr;

	static int ers[8192][8192], ecs[8192][8192], ens[8192];

	int r1, c1, r2, c2;

	//
	for(k=0; k<n; ++k)
	{
		ens[k] = 0;

		//
		for(i=rs[inds[k]]-ss[inds[k]]/2; i<rs[inds[k]]+ss[inds[k]]/2; ++i)
			for(j=cs[inds[k]]-ss[inds[k]]/2; j<cs[inds[k]]+ss[inds[k]]/2; ++j)
				if(edgemaps[inds[k]][i*ldims[inds[k]]+j])
				{
					if( (i-rs[inds[k]])*(i-rs[inds[k]]) + (j-cs[inds[k]])*(j-cs[inds[k]]) < (ss[inds[k]]-2)*(ss[inds[k]]-2)/4 )
					{
						ers[k][ens[k]] = i;
						ecs[k][ens[k]] = j;

						++ens[k];
					}
				}
	}

	//
	newstacksize = stacksize;

	stackbyteptr = (int8_t*)&stack[0];

	if(1)
	{
		maxnumiters = 2048;

		numiters = 0;

		while(newstacksize<maxstacksize && numiters<maxnumiters)
		{
			float o;

			int e;

			//
			k = mwcrand()%n;
			e = mwcrand()%ens[k];

			//
			p = (int)( ss[inds[k]]*s2p );

			//
			o = getorient(ers[k][e], ecs[k][e], pixelss[inds[k]], nrowss[inds[k]], ncolss[inds[k]], ldims[inds[k]], 3);

			//
			r1 = MIN(MAX(rs[inds[k]]-ss[inds[k]]/2+1, ers[k][e]-sin(o)*p), rs[inds[k]]+ss[inds[k]]/2-1);
			c1 = MIN(MAX(cs[inds[k]]-ss[inds[k]]/2+1, ecs[k][e]-cos(o)*p), cs[inds[k]]+ss[inds[k]]/2-1);

			r2 = MIN(MAX(rs[inds[k]]-ss[inds[k]]/2+1, ers[k][e]+sin(o)*p), rs[inds[k]]+ss[inds[k]]/2-1);
			c2 = MIN(MAX(cs[inds[k]]-ss[inds[k]]/2+1, ecs[k][e]+cos(o)*p), cs[inds[k]]+ss[inds[k]]/2-1);

			//
			stackbyteptr[4*newstacksize+0] = _FIXED_POINT_SCALE_*(r1-rs[inds[k]])/ss[inds[k]];
			stackbyteptr[4*newstacksize+1] = _FIXED_POINT_SCALE_*(c1-cs[inds[k]])/ss[inds[k]];
			stackbyteptr[4*newstacksize+2] = _FIXED_POINT_SCALE_*(r2-rs[inds[k]])/ss[inds[k]];
			stackbyteptr[4*newstacksize+3] = _FIXED_POINT_SCALE_*(c2-cs[inds[k]])/ss[inds[k]];

			//
			int ok = 1;

			for(i=0; i<newstacksize; ++i)
				/*
					proximity requirements
				*/
				if(get_bintest_proximity(stack[i], stack[newstacksize]) < (maxnumiters-numiters)*255/maxnumiters)
					ok = 0;

			int nfails = 0;

			for(k=0; k<n; ++k)
			{
				/*
					stability requirements
				*/

				int T[6], Tp[6];

				// compute the transformation matrix
				T[0] = _FIXED_POINT_SCALE_*ss[k]; T[1] = 0; T[2] = _SQR_FIXED_POINT_SCALE_*rs[k];
				T[3] = 0; T[4] = _FIXED_POINT_SCALE_*ss[k]; T[5] = _SQR_FIXED_POINT_SCALE_*cs[k];

				//
				if(0==bintest(stack[newstacksize], threshold, T, pixelss[inds[k]], nrowss[inds[k]], ncolss[inds[k]], ldims[inds[k]]))
					nfails += 1000;

				for(i=0; i<32; ++i)
				{
					Tp[0] = _FIXED_POINT_SCALE_*ss[k]; Tp[1] = 0; Tp[2] = _SQR_FIXED_POINT_SCALE_*(rs[inds[k]]+mwcrand()%(p/2+1)-(p/2));
					Tp[3] = 0; Tp[4] = _FIXED_POINT_SCALE_*ss[k]; Tp[5] = _SQR_FIXED_POINT_SCALE_*(cs[inds[k]]+mwcrand()%(p/2+1)-(p/2));

					if( 0==bintest(stack[newstacksize], threshold, Tp, pixelss[inds[k]], nrowss[inds[k]], ncolss[inds[k]], ldims[inds[k]]) )
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
	}

	//
	return newstacksize;
}

float get_similarity
	(
		int32_t t1[], int r1, int c1, int s1, uint8_t p1[], int nrows1, int ncols1, int ldim1,
		int32_t t2[], int r2, int c2, int s2, uint8_t p2[], int nrows2, int ncols2, int ldim2
	)
{
	int n1, n2;
	int s12, s21;

	//
	n1 = t1[0];
	n2 = t2[0];

	//
	match_template_at(t1, THRESHOLD, r2, c2, s2, &s12, n1, n1, p2, nrows2, ncols2, ldim2);
	match_template_at(t2, THRESHOLD, r1, c1, s1, &s21, n2, n2, p1, nrows1, ncols1, ldim1);

	//
	return ( s12/(float)n1 + s21/(float)n2 )/2.0f;
}

int partition_data(int32_t* templates[], int rs[], int cs[], int ss[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n)
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
				templates[inds[0]], rs[inds[0]], cs[inds[0]], ss[inds[0]], pixelss[inds[0]], nrowss[inds[0]], ncolss[inds[0]], ldims[inds[0]],
				templates[inds[i]], rs[inds[i]], cs[inds[i]], ss[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]
			)
				>=
			get_similarity
			(
				templates[inds[i]], rs[inds[i]], cs[inds[i]], ss[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]],
				templates[inds[n-1]], rs[inds[n-1]], cs[inds[n-1]], ss[inds[n-1]], pixelss[inds[n-1]], nrowss[inds[n-1]], ncolss[inds[n-1]], ldims[inds[n-1]]
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
				templates[inds[0]], rs[inds[0]], cs[inds[0]], ss[inds[0]], pixelss[inds[0]], nrowss[inds[0]], ncolss[inds[0]], ldims[inds[0]],
				templates[inds[j]], rs[inds[j]], cs[inds[j]], ss[inds[j]], pixelss[inds[j]], nrowss[inds[j]], ncolss[inds[j]], ldims[inds[j]]
			)
				<=
			get_similarity
			(
				templates[inds[j]], rs[inds[j]], cs[inds[j]], ss[inds[j]], pixelss[inds[j]], nrowss[inds[j]], ncolss[inds[j]], ldims[inds[j]],
				templates[inds[n-1]], rs[inds[n-1]], cs[inds[n-1]], ss[inds[n-1]], pixelss[inds[n-1]], nrowss[inds[n-1]], ncolss[inds[n-1]], ldims[inds[n-1]]
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

tnode* grow_subtree(int depth, int32_t stack[], int stacksize, int maxnumtests, int rs[], int cs[], int ss[], int32_t* templates[], uint8_t* pixelss[], uint8_t* edgemaps[], int nrowss[], int ncolss[], int ldims[], int inds[], int n)
{
	int i, newstacksize, n1, n2;

	tnode* root = 0;

	//
	if(n==0)
		return 0;

	root = (tnode*)malloc(sizeof(tnode));

	//
	newstacksize = learn_cluster_features(stack, stacksize, stacksize+maxnumtests, 1.5f*S2P, rs, cs, ss, pixelss, edgemaps, nrowss, ncolss, ldims, inds, n, THRESHOLD);

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

	n1 = partition_data(templates, rs, cs, ss, pixelss, nrowss, ncolss, ldims, inds, n); // 0th and (n-1)st samples serve as "anchors"

	n1 = MAX(n1, 1); // hack?

	n2 = n - n1;

	//
	root->subtree1 = grow_subtree(depth+1, stack, newstacksize, maxnumtests, rs, cs, ss, templates, pixelss, edgemaps, nrowss, ncolss, ldims, &inds[0 ], n1);
	root->subtree2 = grow_subtree(depth+1, stack, newstacksize, maxnumtests, rs, cs, ss, templates, pixelss, edgemaps, nrowss, ncolss, ldims, &inds[n1], n2);

	//
	return root;
}

tnode* grow_tree(int rs[], int cs[], int ss[], uint8_t* pixelss[], uint8_t* edgess[], int nrowss[], int ncolss[], int ldims[], int n)
{
	int i;
	int* inds;

	int32_t** templates;

	int maxstacksize;
	int32_t* stack;

	tnode* root;

	//
	inds = (int*)malloc(n*sizeof(int));

	for(i=0; i<n; ++i)
		inds[i] = i;

	// learn templates ... (just for similarity estimation)
	templates = (int32_t**)malloc(n*sizeof(int32_t*));

	for(i=0; i<n; ++i)
	{
		templates[i] = (int32_t*)malloc((MAXNUMTESTS+1)*sizeof(int32_t));
		learn_template(templates[i], MAXNUMTESTS, 1, S2P, rs[i], cs[i], ss[i], pixelss[i], edgess[i], nrowss[i], ncolss[i], ncolss[i], THRESHOLD);
	}

	//
	maxstacksize = n*MAXNUMTESTS;

	stack = (int32_t*)malloc(maxstacksize*sizeof(int32_t));

	//
	root = grow_subtree(0, stack, 0, MAXNUMTESTS/4, rs, cs, ss, templates, pixelss, edgess, nrowss, ncolss, ldims, inds, n);

	//
	for(i=0; i<n; ++i)
		free(templates[i]);
	free(templates);

	free(stack);

	//
	return root;
}

int numtags;
int tags[8192];

int get_tree_output(tnode* root, int threshold, int n0max, int r, int c, int s, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	if(root->template)
	{
		int n1;

		if(!match_template_at(root->template, threshold, r, c, s, &n1, n0max, n0max, pixels, nrows, ncols, ldim))
			return 0;
	}

	if(root->leaf)
	{
		tags[numtags] = root->tag;
		++numtags;

		return 1;
	}
	else
	{
		return
			get_tree_output(root->subtree1, threshold, n0max, r, c, s, pixels, nrows, ncols, ldim)
				|
			get_tree_output(root->subtree2, threshold, n0max, r, c, s, pixels, nrows, ncols, ldim);
	}
}