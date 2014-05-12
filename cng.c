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

int learn_cluster_features(int32_t stack[], int stacksize, int maxstacksize, int* Ts[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n, int32_t tcodepool[], int tcodepoolsize)
{
	int i, j, k, newstacksize, numiters, maxnumiters;

	//
	newstacksize = stacksize;

	if(1)
	{
		maxnumiters = 2048;

		numiters = 0;

		while(newstacksize<maxstacksize && numiters<maxnumiters)
		{
			//
			k = mwcrand()%tcodepoolsize;

			stack[newstacksize] = tcodepool[k];

			//
			/*
			o = get_area_orientation(ers[k][e], ecs[k][e], pixelss[inds[k]], nrowss[inds[k]], ncolss[inds[k]], ldims[inds[k]], 3);

			r1 = ers[k][e]-sin(o)*p;
			c1 = ecs[k][e]-cos(o)*p;

			r2 = ers[k][e]+sin(o)*p;
			c2 = ecs[k][e]+cos(o)*p;

			int det = Ts[k][0]*Ts[k][4] - Ts[k][1]*Ts[k][3];

			stackbyteptr[4*newstacksize+0] = ( +(_SQR_FIXED_POINT_SCALE_*r1-Ts[k][2])*Ts[k][4] + -(_SQR_FIXED_POINT_SCALE_*c1-Ts[k][5])*Ts[k][1] )/det;
			stackbyteptr[4*newstacksize+1] = ( -(_SQR_FIXED_POINT_SCALE_*r1-Ts[k][2])*Ts[k][3] + +(_SQR_FIXED_POINT_SCALE_*c1-Ts[k][5])*Ts[k][0] )/det;
			stackbyteptr[4*newstacksize+2] = ( +(_SQR_FIXED_POINT_SCALE_*r2-Ts[k][2])*Ts[k][4] + -(_SQR_FIXED_POINT_SCALE_*c2-Ts[k][5])*Ts[k][1] )/det;
			stackbyteptr[4*newstacksize+3] = ( -(_SQR_FIXED_POINT_SCALE_*r2-Ts[k][2])*Ts[k][3] + +(_SQR_FIXED_POINT_SCALE_*c2-Ts[k][5])*Ts[k][0] )/det;
			*/

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

				int* T = Ts[k];

				//
				if(0==bintest(stack[newstacksize], THRESHOLD, T, pixelss[inds[k]], nrowss[inds[k]], ncolss[inds[k]], ldims[inds[k]]))
					nfails += 1000;

				for(i=0; i<32; ++i)
				{
					int Tp[6];

					Tp[0] = T[0]; Tp[1] = T[1]; Tp[2] = T[2];
					Tp[3] = T[3]; Tp[4] = T[4]; Tp[5] = T[5];

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

	// compute the transformation matrices
	///match_template_at(t1, THRESHOLD, T2, &s12, n1, n1, p2, nrows2, ncols2, ldim2);
	///match_template_at(t2, THRESHOLD, T1, &s21, n2, n2, p1, nrows1, ncols1, ldim1);

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

tnode* grow_subtree(int depth, int32_t stack[], int stacksize, int maxnumtests, int* Ts[], int32_t* templates[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n, int32_t tcodepool[], int tcodepoolsize)
{
	int i, newstacksize, n1, n2;

	tnode* root = 0;

	//
	if(n==0)
		return 0;

	root = (tnode*)malloc(sizeof(tnode));

	//
	newstacksize = learn_cluster_features(stack, stacksize, stacksize+maxnumtests, Ts, pixelss, nrowss, ncolss, ldims, inds, n, tcodepool, tcodepoolsize);

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
	root->subtree1 = grow_subtree(depth+1, stack, newstacksize, maxnumtests, Ts, templates, pixelss, nrowss, ncolss, ldims, &inds[0 ], n1, tcodepool, tcodepoolsize);
	root->subtree2 = grow_subtree(depth+1, stack, newstacksize, maxnumtests, Ts, templates, pixelss, nrowss, ncolss, ldims, &inds[n1], n2, tcodepool, tcodepoolsize);

	//
	return root;
}

tnode* grow_tree(int* Ts[], int32_t* templates[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int n, int32_t tcodepool[], int tcodepoolsize)
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
	root = grow_subtree(0, stack, 0, MAXNUMTESTS/4, Ts, templates, pixelss, nrowss, ncolss, ldims, inds, n, tcodepool, tcodepoolsize);

	//
	free(stack);

	//
	return root;
}

int numtags;
int tags[8192];

int get_tree_output(tnode* root, int threshold, int n0max, int* T, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	if(root->template)
	{
		int n1;

		if(!match_template_at(root->template, threshold, T, &n1, n0max, pixels, nrows, ncols, ldim))
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
			get_tree_output(root->subtree1, threshold, n0max, T, pixels, nrows, ncols, ldim)
				|
			get_tree_output(root->subtree2, threshold, n0max, T, pixels, nrows, ncols, ldim);
	}
}