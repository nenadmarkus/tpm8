#include <stdio.h>

#include <cv.h>
#include <highgui.h>

/*
	prameters ...
*/

#define MAXNUMTESTS (512)

#define S2P (1.0f/20.0f)
#define THRESHOLD 20

int n0max = 2;
int r0max = 3;

// -----------------------

#define MAXNUMTEMPLATES 8192
int numtemplates = 0;
int32_t templates[MAXNUMTEMPLATES][MAXNUMTESTS+1];
int32_t smoothnesstemplates[MAXNUMTEMPLATES][MAXNUMTESTS+1];

int numtemplateclusters = 0;
int32_t clustertemplates[MAXNUMTEMPLATES][1+MAXNUMTESTS];

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

tnode* root = 0;

float get_overlap(int r1, int c1, int s1, int r2, int c2, int s2)
{
	float overr, overc;

	//
	overr = MAX(0, MIN(r1+s1/2, r2+s2/2) - MAX(r1-s1/2, r2-s2/2));
	overc = MAX(0, MIN(c1+s1/2, c2+s2/2) - MAX(c1-s1/2, c2-s2/2));

	//
	return overr*overc/(float)(s1*s1+s2*s2-overr*overc);
}

void ccdfs(int a[], int i, int rs[], int cs[], int ss[], int n)
{
	int j;

	//
	for(j=0; j<n; ++j)
		if(a[j]==0 && get_overlap(rs[i], cs[i], ss[i], rs[j], cs[j], ss[j])>0.75f)
		{
			//
			a[j] = a[i];

			//
			ccdfs(a, j, rs, cs, ss, n);
		}
}

int find_connected_components(int a[], int rs[], int cs[], int ss[], int n)
{
	int i, ncc, cc;

	//
	if(!n)
		return 0;

	//
	for(i=0; i<n; ++i)
		a[i] = 0;

	//
	ncc = 0;
	cc = 1;

	for(i=0; i<n; ++i)
		if(a[i] == 0)
		{
			//
			a[i] = cc;

			//
			ccdfs(a, i, rs, cs, ss, n);

			//
			++ncc;
			++cc;
		}

	//
	return ncc;
}

int cluster_detections(int a[], int rs[], int cs[], int ss[], int qs[], int n)
{
	int ncc, cc;

	//
	ncc = find_connected_components(a, rs, cs, ss, n);

	if(!ncc)
		return 0;

	//
	int idx = 0;

	for(cc=1; cc<=ncc; ++cc)
	{
		int i, k;

		int sumrs=0, sumcs=0, sumss=0, sumqs=0;

		//
		k = 0;

		for(i=0; i<n; ++i)
			if(a[i] == cc)
			{
				sumrs += rs[i];
				sumcs += cs[i];
				sumss += ss[i];
				sumqs += qs[i];

				++k;
			}
		
		///for(i=0; i<n; ++i)
		///	if(a[i] == cc)
		///		qs[i] = sumqs;

		//
		qs[idx] = sumqs; // accumulated confidence measure

		//
		rs[idx] = sumrs/k;
		cs[idx] = sumcs/k;
		ss[idx] = sumss/k;

		//
		++idx;
	}

	//
	return idx;
}

int match_templates(int rs[], int cs[], int ss[], int qs[], int32_t* ptrs[], int maxndetections,
					uint8_t pixels[], int nrows, int ncols, int ldim,
					float scalefactor, float stridefactor, float minsize, float maxsize, int n0max, int r0max)
{
	int s;
	int ndetections;

	int k = 0;

	int numtags;
	int tags[8192];

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
				int lutidx, i, n1, pass, T[6];

				//
				compute_rcs_transformation(T, r, c, s);

				//
				//*
				numtags = 0;
				get_tree_output(root, THRESHOLD, 3, T, pixels, nrows, ncols, ldim, tags, &numtags);

				for(i=0; i<numtags; ++i)
				{
					pass = match_template_at(templates[tags[i]], THRESHOLD, T, &n1, n0max, pixels, nrows, ncols, ldim);

					if(pass)
					{
						int _n1;

						match_template_at(smoothnesstemplates[tags[i]], THRESHOLD, T, &_n1, MAXNUMTESTS, pixels, nrows, ncols, ldim);

						if(_n1 > smoothnesstemplates[tags[i]][0]/2)
							pass = 0;
					}

					if(pass)
					{
						if(ndetections < maxndetections)
						{
							rs[ndetections] = r;
							cs[ndetections] = c;
							ss[ndetections] = s;

							qs[ndetections] = n1;

							ptrs[ndetections] = templates[tags[i]];

							++ndetections;
						}
					}
				}
				//*/

				//
				/*
				for(i=0; i<numtemplates; ++i)
				{
					pass = match_template_at(templates[i], THRESHOLD, T, &n1, n0max, pixels, nrows, ncols, ldim);

					if(pass)
					{
						int _n1;

						match_template_at(smoothnesstemplates[i], THRESHOLD, T, &_n1, MAXNUMTESTS, pixels, nrows, ncols, ldim);

						if(_n1 > smoothnesstemplates[tags[i]][0]/2)
							pass = 0;
					}

					if(pass)
					{
						if(ndetections < maxndetections)
						{
							rs[ndetections] = r;
							cs[ndetections] = c;
							ss[ndetections] = s;

							qs[ndetections] = n1;

							ptrs[ndetections] = templates[i];

							++ndetections;
						}
					}
				}
				//*/
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
		r1 = (_FIXED_POINT_SCALE_*r + ptr[4*i+0]*s)/_FIXED_POINT_SCALE_;
		c1 = (_FIXED_POINT_SCALE_*c + ptr[4*i+1]*s)/_FIXED_POINT_SCALE_;

		r2 = (_FIXED_POINT_SCALE_*r + ptr[4*i+2]*s)/_FIXED_POINT_SCALE_;
		c2 = (_FIXED_POINT_SCALE_*c + ptr[4*i+3]*s)/_FIXED_POINT_SCALE_;

		//
		cvCircle(drawto, cvPoint((c1+c2)/2, (r1+r2)/2), 1, CV_RGB(0, 255, 0), 2, 8, 0);

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
	int a[MAXNDETECTIONS], rs[MAXNDETECTIONS], cs[MAXNDETECTIONS], ss[MAXNDETECTIONS], qs[MAXNDETECTIONS];
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

	//
	float SCALEFACTOR = 1.1f;
	float STRIDEFACTOR = 0.06f;

	int MINSIZE = 100;//75;
	int MAXSIZE = 600;//150;

	t = getticks();
	ndetections = match_templates(rs, cs, ss, qs, ptrs, MAXNDETECTIONS, pixels, nrows, ncols, ldim, SCALEFACTOR, STRIDEFACTOR, MINSIZE, MAXSIZE, n0max, r0max);

//#define USE_CLUSTERING
#ifdef USE_CLUSTERING
	ndetections = cluster_detections(a, rs, cs, ss, qs, ndetections);
#endif

	t = getticks()-t;

	//
	if(draw)
		for(i=0; i<ndetections; ++i)
		{
#ifndef USE_CLUSTERING
			draw_template_pattern(img, ptrs[i], rs[i], cs[i], ss[i], pixels, nrows, ncols, ldim);
#endif
			cvCircle(img, cvPoint(cs[i], rs[i]), ss[i]/2, CV_RGB(255, 0, 0), 2, 8, 0);
		}

	// if the flag is set, print the results to standard output
	if(print)
	{
		printf("%f [ms] ...\n", 1000.0f*t);

		//for(i=0; i<ndetections; ++i)
		//	printf("%d %d %d %d\n", rs[i], cs[i], ss[i], qs[i]);
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
	int i;

	//
	if(argc >= 2)
	{
		FILE* file = fopen(argv[1], "rb");

		if(!file)
			return 0;

		//
		fread(&numtemplates, sizeof(int), 1, file);

		for(i=0; i<numtemplates; ++i)
		{
			LOAD_TEMPLATE(templates[i], file);
			LOAD_TEMPLATE(smoothnesstemplates[i], file);
		}

		//printf("%d templates loaded ...\n", numtemplates);

		//
		fread(&numtemplateclusters, sizeof(int), 1, file);

		for(i=0; i<numtemplateclusters; ++i)
		{
			LOAD_TEMPLATE(clustertemplates[i], file);
		}

		//
		root = load_tree_from_file(file);

		//
		fclose(file);
	}
	else
		return 0;

	//
	if(argc==3 || argc==4)
	{
		IplImage* img = cvLoadImage(argv[2], CV_LOAD_IMAGE_COLOR);

		process_image(img, 1, 1);

		if(argc==3)
		{
			//*
			cvShowImage("rnt", img);
			cvWaitKey(0);
			//*/
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
