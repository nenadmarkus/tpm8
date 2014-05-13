#include <stdio.h>

#include <cv.h>
#include <highgui.h>

/*
	parameters ...
*/

#define THRESHOLD 25
#define MAXNUMTESTS (128)
#define S2P (1/20.0f)

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

#include "bnt.c"
#include "tme.c"
#include "cng.c"

void draw_template_pattern(IplImage* drawto, int32_t template[], int r, int c, int s, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int i;

	int8_t* ptr;

	//
	ptr = (int8_t*)&template[1];

	for(i=0; i<template[0]; ++i)
	{
		int r1, c1, r2, c2, t;

		//
		r1 = (_FIXED_POINT_SCALE_*r + ptr[4*i+0]*s)/_FIXED_POINT_SCALE_;
		c1 = (_FIXED_POINT_SCALE_*c + ptr[4*i+1]*s)/_FIXED_POINT_SCALE_;

		r2 = (_FIXED_POINT_SCALE_*r + ptr[4*i+2]*s)/_FIXED_POINT_SCALE_;
		c2 = (_FIXED_POINT_SCALE_*c + ptr[4*i+3]*s)/_FIXED_POINT_SCALE_;

		//
		cvLine(drawto, cvPoint(c1, r1), cvPoint(c2, r2), CV_RGB(255, 255, 255), 0, 8, 0);
	}
}

/*
	
*/

#define MAXNUMTEMPLATES 8192
int numtemplates = 0;
int32_t templates[MAXNUMTEMPLATES][MAXNUMTESTS+1];
int32_t smoothnesstemplates[MAXNUMTEMPLATES][MAXNUMTESTS+1];

int numtemplateclusters = 0;
int32_t clustertemplates[MAXNUMTEMPLATES][1+MAXNUMTESTS];

tnode* root = 0;

void learn_templates(uint8_t* pix[], int rs[], int cs[], int ss[], int nrowss[], int ncolss[], int ldims[], int numsamples)
{
	int i, j, n;

	float t;

	uint8_t* edgess[2048];

	int32_t* ptemplates[8192];

	int* Ts[8192];

	//
	t = getticks();

	numtemplates = 0;

	for(n=0; n<numsamples; ++n)
	{
		int lutidx, learnnew;

		//
		learnnew = 1;

		if(learnnew && numtemplates<MAXNUMTEMPLATES)
		{
			IplImage* edges;
			IplImage* img;

			uint8_t* edgemap;

			//
			edgemap = (uint8_t*)malloc( nrowss[n]*ldims[n]*sizeof(uint8_t) );

			//
			img = cvCreateImageHeader(cvSize(ncolss[n], nrowss[n]), IPL_DEPTH_8U, 1);
			img->imageData = (char*)pix[n];
			img->widthStep = ldims[n];

			edges = cvCreateImageHeader(cvSize(ncolss[n], nrowss[n]), IPL_DEPTH_8U, 1);
			edges->imageData = (char*)edgemap;
			edges->widthStep = ldims[n];

			//
			cvCanny(img, edges, 150, 225, 3);

			//cvShowImage("...", edges); cvWaitKey(0);

			//
			learn_template(templates[numtemplates], MAXNUMTESTS, 1, S2P, rs[n], cs[n], ss[n], pix[n], edgemap, nrowss[n], ncolss[n], ldims[n], THRESHOLD);
			learn_template(smoothnesstemplates[numtemplates], MAXNUMTESTS, 0, S2P, rs[n], cs[n], ss[n], pix[n], pix[n], nrowss[n], ncolss[n], ldims[n], THRESHOLD);

			///draw_template_pattern(edges, templates[numtemplates], rs[n], cs[n], ss[n], pix[n], nrowss[n], ncolss[n], ldims[n]); cvCircle(edges, cvPoint(cs[n], rs[n]), ss[n]/2, CV_RGB(255, 255, 255), 2, 8, 0); cvShowImage("...", edges); cvWaitKey(0);

			edgess[numtemplates] = edgemap;

			//
			Ts[numtemplates] = (int*)malloc(6*sizeof(int));
			compute_rcs_transformation(Ts[numtemplates], rs[n], cs[n], ss[n]);

			ptemplates[numtemplates] = &templates[numtemplates][0];

			//
			++numtemplates;

			//
			///free(edgemap);

			cvReleaseImageHeader(&edges);
			cvReleaseImageHeader(&img);
		}
	}

	printf("%d templates learned in %f [ms]\n", numtemplates, 1000.0f*(getticks()-t));

	//
	int tcodepoolsize;
	int32_t tcodepool[1024*MAXNUMTESTS], tmptemplate[MAXNUMTESTS+1];

	tcodepoolsize = 0;

	for(i=0; i<numtemplates; ++i)
	{
		learn_template(tmptemplate, MAXNUMTESTS, 1, 1.5f*S2P, rs[i], cs[i], ss[i], pix[i], edgess[i], nrowss[i], ncolss[i], ldims[i], THRESHOLD);

		for(j=0; j<tmptemplate[0]; ++j)
		{
			tcodepool[tcodepoolsize] = tmptemplate[j+1];

			++tcodepoolsize;
		}
	}

	int perturbationstrength = (int)( 1.5f*S2P*ss[0]/2 );

	t = getticks();
	root = grow_tree(Ts, ptemplates, pix, nrowss, ncolss, ldims, numsamples, tcodepool, tcodepoolsize, perturbationstrength);
	printf("%f [ms] elapsed for clustering\n", 1000.0f*(getticks()-t));
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

int load_samples(char* folder, uint8_t* pix[], int rs[], int cs[], int ss[], int nrowss[], int ncolss[], int ldims[], int maxn)
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
		IplImage* tmp;
		uint8_t* p = 0;

		sprintf(path, "%s/%s", folder, name);

		tmp = cvLoadImage(path, CV_LOAD_IMAGE_GRAYSCALE);

		///if( loadrid(&p, &nrowss[n], &ncolss[n], path) )
		if(tmp)
		{
			//
			pix[n] = (uint8_t*)tmp->imageData;

			nrowss[n] = tmp->height;
			ncolss[n] = tmp->width;
			ldims[n] = tmp->widthStep;

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
	int n, tdepth, i;

	#define MAXN 1024

	static uint8_t* pix[MAXN];
	static int rs[MAXN], cs[MAXN], ss[MAXN], nrowss[MAXN], ncolss[MAXN], ldims[MAXN];

	//
	smwcrand(time(0));
	///smwcrand(12344);

	//
	if(argc != 3)
		return 0;

	n = load_samples(argv[1], pix, rs, cs, ss, nrowss, ncolss, ldims, MAXN);

	//
	t = getticks();
	learn_templates(pix, rs, cs, ss, nrowss, ncolss, ldims, n);
	printf("elapsed time: %f [s]\n", getticks()-t);

	//
	if(numtemplates)
	{
		FILE* file = fopen(argv[2], "wb");

		if(!file)
		{
			printf("cannot save results to '%s'", argv[3]);
			return 0;
		}

		//
		fwrite(&numtemplates, sizeof(int), 1, file);

		for(i=0; i<numtemplates; ++i)
		{
			SAVE_TEMPLATE(templates[i], file);
			SAVE_TEMPLATE(smoothnesstemplates[i], file);
		}

		//
		fwrite(&numtemplateclusters, sizeof(int), 1, file);

		for(i=0; i<numtemplateclusters; ++i)
		{
			SAVE_TEMPLATE(clustertemplates[i], file);
		}

		//
		save_tree_to_file(root, file);

		//
		fclose(file);
	}

	//
	return 0;
}