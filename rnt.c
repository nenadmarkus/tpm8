#include <stdio.h>

#include <cv.h>
#include <highgui.h>

/*
	prameters ...
*/

#define USE_CLAHE

#define NTESTS (256)

#define THRESHOLD 25

int n0max = 3;
int r0max = 3;

// -----------------------

#define TDEPTH 12

int32_t tree[1<<(TDEPTH+1)];

int templatecounts[1<<TDEPTH];
int32_t templatelut[1<<TDEPTH][64][1+NTESTS];

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

#define USE_RGB

#include "bnt.c"
#include "tme.c"
#include "cng.c"
#include "ahe.c"

int match_templates(int rs[], int cs[], int ss[], int32_t* ptrs[], int maxndetections,
					uint8_t pixels[], int nrows, int ncols, int ldim,
					float scalefactor, float stridefactor, float minsize, float maxsize, int n0max, int r0max)
{
	int s;
	int ndetections;

	//
	ndetections = 0;

	s = minsize;

	while(s<=maxsize)
	{
		int r, c, dr, dc;

		//
		dr = MAX(stridefactor*s, 1);
		dc = MAX(stridefactor*s, 1);

		//
		for(r=s/2+1; r<=nrows-s/2-1; r+=dr)
			for(c=s/2+1; c<=ncols-s/2-1; c+=dc)
			{
				int lutidx, i, n1, pass;

				//
				lutidx = get_tree_output(tree, THRESHOLD, r, c, s, pixels, nrows, ncols, ldim);

				//
				for(i=0; i<templatecounts[lutidx]; ++i)
				{
					int32_t* template = (int32_t*)&templatelut[lutidx][i][0];

					pass = match_template_at(template, THRESHOLD, r, c, s, &n1, n0max, r0max, pixels, nrows, ncols, ldim);

					if(pass)
					{
						if(ndetections < maxndetections)
						{
							rs[ndetections] = r;
							cs[ndetections] = c;
							ss[ndetections] = s;

							ptrs[ndetections] = templatelut[lutidx][i];

							++ndetections;
						}
					}
				}
			}

		//
		s = MAX(scalefactor*s, s+1);
	}

	//
	return ndetections;
}

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
		r1 = (NORMALIZATION*r + ptr[4*i+0]*s)/NORMALIZATION;
		c1 = (NORMALIZATION*c + ptr[4*i+1]*s)/NORMALIZATION;

		r2 = (NORMALIZATION*r + ptr[4*i+2]*s)/NORMALIZATION;
		c2 = (NORMALIZATION*c + ptr[4*i+3]*s)/NORMALIZATION;

		//
		//*
		if(ABS(pixels[r1*ldim+c1]-pixels[r2*ldim+c2]) > THRESHOLD)
			cvCircle(drawto, cvPoint((c1+c2)/2, (r1+r2)/2), 1, CV_RGB(0, 255, 0), 2, 8, 0);
		else
			cvCircle(drawto, cvPoint((c1+c2)/2, (r1+r2)/2), 1, CV_RGB(255, 0, 0), 2, 8, 0);
		//*/

		/*
		if(ABS(pixels[r1*ldim+c1]-pixels[r2*ldim+c2]) > THRESHOLD)
			cvLine(drawto, cvPoint(c1, r1), cvPoint(c2, r2), CV_RGB(0, 255, 0), 0, 8, 0);
		else
			cvLine(drawto, cvPoint(c1, r1), cvPoint(c2, r2), CV_RGB(255, 0, 0), 0, 8, 0);
		//*/
	}
}

/*
	
*/

void process_image(IplImage* img, int draw, int print)
{
	int i, j;

	uint8_t* pixels;
	int nrows, ncols, ldim;

	static IplImage* toprocess = 0;

	#define MAXNDETECTIONS 8192
	int ndetections;
	int rs[MAXNDETECTIONS], cs[MAXNDETECTIONS], ss[MAXNDETECTIONS];
	int32_t* ptrs[MAXNDETECTIONS];

	float t;

#ifndef USE_RGB
	// grayscale image
	if(!toprocess)
		toprocess = cvCreateImage(cvSize(img->width, img->height), img->depth, 1);
	if(img->nChannels == 3)
		cvCvtColor(img, toprocess, CV_RGB2GRAY);
	else
		cvCopy(img, toprocess, 0);
#else
	toprocess = img;
#endif

	// get relevant image data
	pixels = (uint8_t*)toprocess->imageData;
	nrows = toprocess->height;
	ncols = toprocess->width;
	ldim = toprocess->widthStep;

#ifndef USE_RGB
#ifdef USE_CLAHE
	CLAHE(pixels, pixels, nrows, ncols, ldim, 8, 8, 3);
#endif
#endif

	//
	float SCALEFACTOR = 1.05f;
	float STRIDEFACTOR = 0.04f;

	int MINSIZE = 70;
	int MAXSIZE = 120;

	t = getticks();
	ndetections = match_templates(rs, cs, ss, ptrs, MAXNDETECTIONS, pixels, nrows, ncols, ldim, SCALEFACTOR, STRIDEFACTOR, MINSIZE, MAXSIZE, n0max, r0max);
	t = getticks()-t;

	//
	if(draw)
		for(i=0; i<ndetections; ++i)
		{
			draw_template_pattern(img, ptrs[i], rs[i], cs[i], ss[i], pixels, nrows, ncols, ldim);

			cvCircle(img, cvPoint(cs[i], rs[i]), ss[i]/2, CV_RGB(255, 0, 0), 2, 8, 0);
		}

	// if the flag is set, print the results to standard output
	if(print)
	{
		printf("%f [ms] ...\n", 1000.0f*t);
	}
}

void process_video_frames(char* path)
{
	const char* windowname = "rnt";

	CvCapture* capture;

	IplImage* frame;
	IplImage* framecopy;

	int stop;

	//
	cvNamedWindow(windowname, 1);

	// try to initialize video capture from the default webcam
	if(!path)
		capture = cvCaptureFromCAM(0);
	else
		capture = cvCaptureFromAVI(path);

	if(!capture)
	{
		printf("Cannot initialize video capture!\n");
		return;
	}

	// start the main loop in which we'll process webcam output
	framecopy = 0;
	stop = 0;
	while(!stop)
	{
		int key;

		// wait ...
		if(!path)
			key = cvWaitKey(5);
		else
			key = cvWaitKey(0);

		// retrieve a pointer to the image acquired from the webcam
		if(!cvGrabFrame(capture))
			break;
		frame = cvRetrieveFrame(capture, 1);

		// we terminate the loop if we don't get any data from the webcam or the user has pressed 'q'
		if(!frame || key=='q')
			stop = 1;
		else if(key=='p')
		{
			n0max += 1;

			printf("---------------------- n0max = %d\n", n0max);
		}
		else if(key=='m')
		{
			n0max -= 1;

			if(n0max < 0)
				n0max = 0;

			printf("---------------------- n0max = %d\n", n0max);
		}
		else
		{
			// we mustn't tamper with internal OpenCV buffers and that's the reason why we're making a copy of the current frame
			if(!framecopy)
				framecopy = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, frame->nChannels);
			cvCopy(frame, framecopy, 0);

			///cvFlip(framecopy, framecopy, 1); // webcam outputs mirrored frames (at least on my machines); you can safely comment out this line if you find it unnecessary

			//
			process_image(framecopy, 1, 1);

			// display the image to the user
			cvShowImage(windowname, framecopy);
		}
	}

	// cleanup
	cvReleaseImage(&framecopy);
	cvReleaseCapture(&capture);
	cvDestroyWindow(windowname);
}

int main(int argc, char* argv[])
{
	int i, j;

	//
	if(argc >= 2)
	{
		FILE* file = fopen(argv[1], "rb");

		if(!file)
			return 0;

		//
		fread(&tree[0], sizeof(int32_t), 1, file);
		fread(&tree[1], sizeof(int32_t), (1<<(tree[0]+1))-1, file);

		//
		for(i=0; i<(1<<tree[0]); ++i)
		{
			fread(&templatecounts[i], sizeof(int32_t), 1, file);

			for(j=0; j<templatecounts[i]; ++j)
			{
				fread(&templatelut[i][j][0], sizeof(int32_t), 1, file);
				fread(&templatelut[i][j][1], sizeof(int32_t), templatelut[i][j][0], file);
			}
		}

		fclose(file);
	}
	else
		return 0;

	printf("tree depth = %d\n", tree[0]); j=0; for(i=0; i<(1<<tree[0]); ++i) j+=templatecounts[i]; printf("template count = %d\n", j);

	//
	if(argc==3 || argc==4)
	{
		IplImage* img = cvLoadImage(argv[2], CV_LOAD_IMAGE_COLOR);

		process_image(img, 1, 1);

		if(argc==3)
		{
			cvShowImage("rnt", img);
			cvWaitKey(0);
		}
		else
			cvSaveImage(argv[3], img, 0);

		return 0;
	}
	else
		process_video_frames(0);

	//
	return 0;
}
