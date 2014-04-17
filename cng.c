
//#define SAVE_TREE(tree, file) fwrite((tree), sizeof(int32_t), 1<<((tree)[0]+1), (file));
//#define LOAD_TREE(tree, file) (tree)=(int32_t*)malloc(sizeof(int32_t)), fread(&(tree)[0], sizeof(int32_t), 1, (file)), (tree)=(int32_t*)realloc((tree), (1<<((tree)[0]+1))*sizeof(int32_t)),fread(&(tree)[1], sizeof(int32_t), 1<<((tree)[0]+1), (file));

#define SAVE_TREE(tree, file) fwrite((tree), sizeof(int32_t), 1<<(tree)[0], (file));
#define LOAD_TREE(tree, file) (tree)=(int32_t*)malloc(sizeof(int32_t)), fread(&(tree)[0], sizeof(int32_t), 1, (file)), (tree)=(int32_t*)realloc((tree), (1<<(tree)[0])*sizeof(int32_t)),fread(&(tree)[1], sizeof(int32_t), (1<<(tree)[0])-1, (file));

float get_test_quality(int tcode, int rs[], int cs[], int ss[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n)
{
	int i, n0, n1;

	//
	n0 = 0;
	n1 = 0;

	for(i=0; i<n; ++i)
		if( !bintest(tcode, THRESHOLD, rs[inds[i]], cs[inds[i]], ss[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
			++n0;
		else
			++n1;

	//
	return MIN(n0, n1);
}

int split_data(int tcode, int rs[], int cs[], int ss[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n)
{
	int stop;
	int i, j;

	int n0;

	//
	stop = 0;

	i = 0;
	j = n - 1;

	while(!stop)
	{
		//
		while( !bintest(tcode, THRESHOLD, rs[inds[i]], cs[inds[i]], ss[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
		{
			if( i==j )
				break;
			else
				++i;
		}

		while( bintest(tcode, THRESHOLD, rs[inds[j]], cs[inds[j]], ss[inds[j]], pixelss[inds[j]], nrowss[inds[j]], ncolss[inds[j]], ldims[inds[j]]) )
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
		{
			// swap
			inds[i] = inds[i] ^ inds[j];
			inds[j] = inds[i] ^ inds[j];
			inds[i] = inds[i] ^ inds[j];
		}
	}

	//
	n0 = 0;

	for(i=0; i<n; ++i)
		if( !bintest(tcode, THRESHOLD, rs[inds[i]], cs[inds[i]], ss[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
			++n0;

	//
	return n0;
}

int grow_subtree(int32_t tcodes[], int nodeidx, int d, int maxd, int rs[], int cs[], int ss[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n)
{
	int i, nrands, n0, n1, Q;

	//
	if(d >= maxd)
	{
		///printf("%d ", n);

		return 1;
	}

	if(n <= 1)
	{
		//
		tcodes[nodeidx] = mwcrand();

		//
		grow_subtree(tcodes, 2*nodeidx+1, d+1, maxd, rs, cs, ss, pixelss, nrowss, ncolss, ldims, inds, n);
		grow_subtree(tcodes, 2*nodeidx+2, d+1, maxd, rs, cs, ss, pixelss, nrowss, ncolss, ldims, inds, n);

		return 1;
	}

	//
	nrands = 128;

	Q = 0;
	tcodes[nodeidx] = mwcrand();

	for(i=0; i<nrands; ++i)
	{
		int q;
		int32_t tcode;

		//
		tcode = mwcrand();

		q = get_test_quality(tcode, rs, cs, ss, pixelss, nrowss, ncolss, ldims, inds, n);

		//
		if(q > Q)
		{
			Q = q;
			tcodes[nodeidx] = tcode;
		}
	}

	//
	n0 = split_data(tcodes[nodeidx], rs, cs, ss, pixelss, nrowss, ncolss, ldims, inds, n);

	n1 = n - n0;

	//
	grow_subtree(tcodes, 2*nodeidx+1, d+1, maxd, rs, cs, ss, pixelss, nrowss, ncolss, ldims, &inds[0 ], n0);
	grow_subtree(tcodes, 2*nodeidx+2, d+1, maxd, rs, cs, ss, pixelss, nrowss, ncolss, ldims, &inds[n0], n1);

	//
	return 1;
}

int32_t* grow_tree(int d, int rs[], int cs[], int ss[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int n)
{
	int i;

	int* inds;

	int32_t* tree;

	//
	inds = (int*)malloc(n*sizeof(int));

	for(i=0; i<n; ++i)
		inds[i] = i;

	//
	///tree = (int32_t*)malloc( (1<<(d+1))*sizeof(int32_t) );
	tree = (int32_t*)malloc( (1<<d)*sizeof(int32_t) );

	//
	tree[0] = d;

	if(!grow_subtree(&tree[1], 0, 0, d, rs, cs, ss, pixelss, nrowss, ncolss, ldims, inds, n))
	{
		free(inds);
		free(tree);

		return 0;
	}
	else
	{
		free(inds);

		return tree;
	}
}

int get_tree_output(int32_t tree[], int threshold, int r, int c, int s, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int d, tdepth, idx;
	int32_t* tcodes;

	//
	tdepth = tree[0];
	tcodes = &tree[1];

	idx = 0;

	for(d=0; d<tdepth; ++d)
	{
		if( bintest(tcodes[idx], threshold, r, c, s, pixels, nrows, ncols, ldim) )
			idx = 2*idx + 2;
		else
			idx = 2*idx + 1;
	}

	//
	return idx - ((1<<tdepth)-1);
}