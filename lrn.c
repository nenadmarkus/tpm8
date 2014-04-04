#include <stdio.h>

#include <cv.h>
#include <highgui.h>

/*
	...
*/

#ifndef MAX
#define MAX(a, b) ((a)>(b)?(a):(b))
#endif

#ifndef MIN
#define MIN(a, b) ((a)<(b)?(a):(b))
#endif

#ifndef ABS
#define ABS(x) ((x)>0?(x):(-(x)))
#endif

/*
	portable time function
*/

#ifdef __GNUC__
#include <time.h>
float getticks()
{
	struct timespec ts;

	if(clock_gettime(CLOCK_MONOTONIC, &ts) < 0)
	{
		printf("clock_gettime error\n");

		return -1.0f;
	}

	return ts.tv_sec + 1e-9f*ts.tv_nsec;
}
#else
#include <windows.h>
float getticks()
{
	static double freq = -1.0;
	LARGE_INTEGER lint;

	if(freq < 0.0)
	{
		if(!QueryPerformanceFrequency(&lint))
			return -1.0f;

		freq = lint.QuadPart;
	}

	if(!QueryPerformanceCounter(&lint))
		return -1.0f;

	return (float)( lint.QuadPart/freq );
}
#endif

/*
	PRNG
*/

uint32_t mwcrand_r(uint64_t* state)
{
	uint32_t* m;

	//
	m = (uint32_t*)state;

	// bad state?
	if(m[0] == 0)
		m[0] = 0xAAAA;

	if(m[1] == 0)
		m[1] = 0xBBBB;

	// mutate state
	m[0] = 36969 * (m[0] & 65535) + (m[0] >> 16);
	m[1] = 18000 * (m[1] & 65535) + (m[1] >> 16);

	// output
	return (m[0] << 16) + m[1];
}

uint64_t prngglobal = 0x12345678000fffffLL;

void smwcrand(uint32_t seed)
{
	prngglobal = 0x12345678000fffffLL*seed;
}

uint32_t mwcrand()
{
	return mwcrand_r(&prngglobal);
}

/*
	
*/

int loadrid(uint8_t* pixels[], int* nrows, int* ncols, const char* path)
{
	FILE* file;
	int w, h;

	// open file
	file = fopen(path, "rb");

	if(!file)
	{
		return 0;
	}

	// read width
	fread(&w, sizeof(int), 1, file);
	// read height
	fread(&h, sizeof(int), 1, file);

	// allocate image memory
	*nrows = h;
	*ncols = w;

	*pixels = (uint8_t*)malloc(w*h*sizeof(uint8_t));

	if(!*pixels)
	{
		fclose(file);

		return 0;
	}

	// read image data
	fread(*pixels, sizeof(uint8_t), w*h, file);

	// clean up
	fclose(file);

	// we're done
	return 1;
}

/*
	
*/

#define NTESTS 96
#define S2P (1/20.0f)

int n0max = 5;

#include "tm.c"

/*
	template clustering
*/

float get_test_quality(int tcode, int rs[], int cs[], int ss[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n)
{
	int i, n0, n1;

	//
	n0 = 0;
	n1 = 0;

	for(i=0; i<n; ++i)
		if( !bintest(tcode, rs[inds[i]], cs[inds[i]], ss[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
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
		while( !bintest(tcode, rs[inds[i]], cs[inds[i]], ss[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
		{
			if( i==j )
				break;
			else
				++i;
		}

		while( bintest(tcode, rs[inds[j]], cs[inds[j]], ss[inds[j]], pixelss[inds[j]], nrowss[inds[j]], ncolss[inds[j]], ldims[inds[j]]) )
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
		if( !bintest(tcode, rs[inds[i]], cs[inds[i]], ss[inds[i]], pixelss[inds[i]], nrowss[inds[i]], ncolss[inds[i]], ldims[inds[i]]) )
			++n0;

	//
	return n0;
}

int grow_subtree(int32_t tcodes[], int nodeidx, int d, int maxd, int rs[], int cs[], int ss[], uint8_t* pixelss[], int nrowss[], int ncolss[], int ldims[], int inds[], int n)
{
	int i, nrands, n0, n1, Q;

	//
	if(d > maxd)
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
	nrands = 1024;

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
	tree = (int32_t*)malloc( (1<<(d+1))*sizeof(int32_t) );

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

/*
	
*/

int32_t* tree = 0;

#define TDEPTH 12

int templatecounts[1<<TDEPTH];
int32_t templatelut[1<<TDEPTH][64][1+NTESTS];

void learn_templates(uint8_t* pix[], int rs[], int cs[], int ss[], int nrowss[], int ncolss[], int ntemplates, int tdepth)
{
	int n;

	//
	tree = grow_tree(tdepth, rs, cs, ss, pix, nrowss, ncolss, ncolss, ntemplates);

	//
	for(n=0; n<(1<<tdepth); ++n) templatecounts[n] = 0;

	for(n=0; n<ntemplates; ++n)
	{
		int i, lutidx, learnnew;

		//
		lutidx = get_tree_output(tree, rs[n], cs[n], ss[n], pix[n], nrowss[n], ncolss[n], ncolss[n]);

		//
		learnnew = 1;

		for(i=0; i<templatecounts[lutidx]; ++i)
		{
			int n1;

			if( match_template_at(templatelut[lutidx][i], rs[n], cs[n], ss[n], &n1, n0max, pix[n], nrowss[n], ncolss[n], ncolss[n]) )
			{
				learnnew = 0;
			}
		}

		if(learnnew)
		{
			IplImage* edges;
			IplImage* img;

			uint8_t* edgemap;

			//
			edgemap = (uint8_t*)malloc( nrowss[n]*ncolss[n]*sizeof(uint8_t) );

			//
			img = cvCreateImageHeader(cvSize(ncolss[n], nrowss[n]), IPL_DEPTH_8U, 1);
			img->imageData = pix[n];
			img->widthStep = ncolss[n];

			edges = cvCreateImageHeader(cvSize(ncolss[n], nrowss[n]), IPL_DEPTH_8U, 1);
			edges->imageData = edgemap;
			edges->widthStep = ncolss[n];

			cvCanny(img, edges, 75, 125, 3);

			//cvShowImage("...", edges); cvWaitKey(0);

			//
			learn_template(templatelut[lutidx][templatecounts[lutidx]], NTESTS, S2P, rs[n], cs[n], ss[n], pix[n], edgemap, nrowss[n], ncolss[n], ncolss[n]);

			++templatecounts[lutidx];

			//
			free(edgemap);

			cvReleaseImage(&edges);
			cvReleaseImageHeader(&img);
		}

		///printf("%d\n", n);
	}
}

/*
	
*/

void display_image(uint8_t pixels[], int nrows, int ncols, int ldim)
{
	IplImage* header;

	//
	header = cvCreateImageHeader(cvSize(ncols, nrows), IPL_DEPTH_8U, 1);

	header->height = nrows;
	header->width = ncols;
	header->widthStep = ldim;

	header->imageData = (char*)pixels;

	//
	cvCircle(header, cvPoint(ncols/2, nrows/2), 2*MIN(nrows, ncols)/3/2, CV_RGB(0, 0, 0), 2, 8, 0);

	//
	cvShowImage("...", header);

	//
	cvReleaseImageHeader(&header);
}

int load_templates(char* folder, uint8_t* pix[], int rs[], int cs[], int ss[], int nrowss[], int ncolss[], int maxn)
{
	int n;

	FILE* list;

	char path[1024], name[1024];

	//
	sprintf(path, "%s/list", folder);

	list = fopen(path, "r");

	if(!list)
	{
		printf("cannot open '%s'", path);

		return 0;
	}

	//
	n = 0;

	while( fscanf(list, "%s", name) == 1 )
	{
		uint8_t* p = 0;

		sprintf(path, "%s/%s", folder, name);

		if( loadrid(&p, &nrowss[n], &ncolss[n], path) )
		{
			//
			pix[n] = p;

			rs[n] = nrowss[n]/2;
			cs[n] = ncolss[n]/2;
			ss[n] = 2*MIN(nrowss[n], ncolss[n])/3;

			///display_image(p, nrowss[n], ncolss[n], ncolss[n]); cvWaitKey(0);

			//
			++n;
		}
	}

	//
	return n;
}

/*
	
*/

int main(int argc, char* argv[])
{
	float t;
	int n, tdepth;

	#define MAXN (1<<16)

	static uint8_t* pix[MAXN];
	static int rs[MAXN], cs[MAXN], ss[MAXN], nrowss[MAXN], ncolss[MAXN];

	//
	if(argc != 4)
		return 0;

	n = load_templates(argv[1], pix, rs, cs, ss, nrowss, ncolss, MAXN);

	sscanf(argv[2], "%d", &tdepth);

	//
	t = getticks();
	learn_templates(pix, rs, cs, ss, nrowss, ncolss, n, tdepth);
	printf("elapsed time: %f [s]\n", getticks()-t);

	//
	if(tree)
	{
		int i, j;

		FILE* file = fopen(argv[3], "wb");

		if(!file)
		{
			printf("cannot save results to '%s'", argv[3]);
			return 0;
		}

		//
		fwrite(tree, sizeof(int32_t), 1<<(tree[0]+1), file);

		//
		for(i=0; i<(1<<tree[0]); ++i)
		{
			fwrite(&templatecounts[i], sizeof(int32_t), 1, file);

			for(j=0; j<templatecounts[i]; ++j)
			{
				fwrite(templatelut[i][j], sizeof(int32_t), 1+templatelut[i][j][0], file);
			}
		}

		fclose(file);

		//
		for(i=0; i<(1<<tree[0]); ++i) printf("%d ", templatecounts[i]); printf("%d\n");
	}

	//
	return 0;
}