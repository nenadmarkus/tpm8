#include <stdio.h>

#include <cv.h>
#include <highgui.h>

/*
	- <RID> file contents:
		- a 32-bit signed integer w (image width)
		- a 32-bit signed integer h (image height)
		- an array of w*h unsigned bytes representing pixel intensities
*/

int saverid(const char* dst, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int r;
	FILE* file;

	//
	file = fopen(dst, "wb");

	if(!file)
		return 0;

	//
	fwrite(&ncols, sizeof(int), 1, file);
	fwrite(&nrows, sizeof(int), 1, file);

	for(r=0; r<nrows; ++r)
		fwrite(&pixels[r*ldim], sizeof(uint8_t), ncols, file);

	//
	fclose(file);

	// we're done
	return 1;
}

/*
	
*/

void warp_template(IplImage* in, int r, int c, int s, IplImage* out)
{
	//
	float tr = s*(0.06f - 0.12f*rand()/(float)RAND_MAX);
	float tc = s*(0.06f - 0.12f*rand()/(float)RAND_MAX);

	float theta = 2*3.14159265*rand()/(float)RAND_MAX;
	float phi = 2*3.14159265*rand()/(float)RAND_MAX;

	float s1 = 0.85f + 0.3f*rand()/(float)RAND_MAX;
	float s2 = 0.85f + 0.3f*rand()/(float)RAND_MAX;

	CvMat cvmat0, cvmat1, cvmat2, cvmat3, cvmat4, cvmat5, H;

	//
	float T0[3][3], T1[3][3], T2[3][3], T3[3][3], T4[3][3], T5[3][3];

	//
	T0[0][0] = 1.0f; T0[0][1] = 0.0f; T0[0][2] = -c + tc;
	T0[1][0] = 0.0f; T0[1][1] = 1.0f; T0[1][2] = -r + tr;
	T0[2][0] = 0.0f; T0[2][1] = 0.0f; T0[2][2] = 1.0f;

	cvmat0 = cvMat(3, 3, CV_32FC1, &T0[0][0]);

	//
	T1[0][0] = cos(phi); T1[0][1] = -sin(phi); T1[0][2] = 0.0f;
	T1[1][0] = sin(phi); T1[1][1] = cos(phi); T1[1][2] = 0.0f;
	T1[2][0] = 0.0f; T1[2][1] = 0.0f; T1[2][2] = 1.0f;

	cvmat1 = cvMat(3, 3, CV_32FC1, &T1[0][0]);

	//
	T2[0][0] = s1; T2[0][1] = 0.0f; T2[0][2] = 0.0f;
	T2[1][0] = 0.0f; T2[1][1] = s2; T2[1][2] = 0.0f;
	T2[2][0] = 0.0f; T2[2][1] = 0.0f; T2[2][2] = 1.0f;

	cvmat2 = cvMat(3, 3, CV_32FC1, &T2[0][0]);

	//
	T3[0][0] = cos(-phi); T3[0][1] = -sin(-phi); T3[0][2] = 0.0f;
	T3[1][0] = sin(-phi); T3[1][1] = cos(-phi); T3[1][2] = 0.0f;
	T3[2][0] = 0.0f; T3[2][1] = 0.0f; T3[2][2] = 1.0f;
	
	cvmat3 = cvMat(3, 3, CV_32FC1, &T3[0][0]);

	//
	T4[0][0] = cos(theta); T4[0][1] = -sin(theta); T4[0][2] = 0.0f;
	T4[1][0] = sin(theta); T4[1][1] = cos(theta); T4[1][2] = 0.0f;
	T4[2][0] = 0.0f; T4[2][1] = 0.0f; T4[2][2] = 1.0f;

	cvmat4 = cvMat(3, 3, CV_32FC1, &T4[0][0]);

	//
	T5[0][0] = 1.0f; T5[0][1] = 0.0f; T5[0][2] = out->width/2;
	T5[1][0] = 0.0f; T5[1][1] = 1.0f; T5[1][2] = out->height/2;
	T5[2][0] = 0.0f; T5[2][1] = 0.0f; T5[2][2] = 1.0f;

	cvmat5 = cvMat(3, 3, CV_32FC1, &T5[0][0]);

	//
	H = cvMat(3, 3, CV_32FC1, &T0[0][0]);

	cvMatMul(&cvmat1, &H, &H);
	cvMatMul(&cvmat2, &H, &H);
	cvMatMul(&cvmat3, &H, &H);
	cvMatMul(&cvmat4, &H, &H);
	cvMatMul(&cvmat5, &H, &H);

	//
	cvZero(out);

	cvWarpPerspective(in, out, &H, 0, cvScalarAll(0));
}

int generate_warps(IplImage* imgs[], int rs[], int cs[], int ss[], int nimgs, int nwarpsperimage, IplImage* warps[])
{
	int i, j, n;

	//
	n = 0;

	for(i=0; i<nimgs; ++i)
	{
		IplImage* template = 0;

		//
		template = cvCreateImage(cvSize(2*ss[i], 2*ss[i]), IPL_DEPTH_8U, 1);

		cvSetImageROI(imgs[i], cvRect(cs[i]-ss[i], rs[i]-ss[i], 2*ss[i], 2*ss[i]));
		cvCopy(imgs[i], template, 0);
		cvResetImageROI(imgs[i]);

		//
		for(j=0; j<nwarpsperimage; ++j)
		{
			int r, c, s;

			//
			r = template->height/2;
			c = template->width/2;
			s = ss[i];

			//
			warps[n] = cvCreateImage(cvSize(template->width, template->height), IPL_DEPTH_8U, 1);

			warp_template(template, r, c, s, warps[n]);

			//
			///cvShowImage("...", warps[n]); cvWaitKey(0);

			//
			++n;
		}

		//
		cvReleaseImage(&template);
	}

	//
	return n;
}

/*
	
*/

const int KEY_SPACE = 32;
const char* windowname = "...";

int mouseeventprocessing = 0;

IplImage* selectin = 0;
IplImage* tmpimg = 0;

int nclicks = 0;
int _r, _c, _s;

void mouse_callback(int e, int x, int y, int flags, void* params)
{
	//
	if(!mouseeventprocessing)
		return;

	if(e==CV_EVENT_LBUTTONDOWN)
	{
		static IplImage* gray = 0;

		//
		if(!gray)
			gray = cvCreateImage(cvSize(selectin->width, selectin->height), selectin->depth, 1);
		if(selectin->nChannels == 3)
			cvCvtColor(selectin, gray, CV_RGB2GRAY);
		else
			cvCopy(selectin, gray, 0);

		//
		if(nclicks==0)
		{
			//
			_r = y;
			_c = x;

			//
			cvCircle(tmpimg, cvPoint(x, y), 1, CV_RGB(255,0,0), 2, 8, 0);

			//
			nclicks = 1;
		}
		else if(nclicks==1)
		{
			//
			_s = 2*sqrt( (_r-y)*(_r-y) + (_c-x)*(_c-x) );

			//
			cvCircle(tmpimg, cvPoint(_c, _r), _s/2, CV_RGB(255,0,0), 2, 8, 0);

			//
			nclicks = 2;
		}

		//
		cvShowImage(windowname, tmpimg);
	}
	else if(e==CV_EVENT_RBUTTONDOWN)
	{
		//
		nclicks = 0;

		//
		cvCopy(selectin, tmpimg, 0);
		cvShowImage(windowname, tmpimg);
	}
}

int select_region(IplImage* frame, int* pr, int* pc, int* ps)
{
	int stop;

	//
	nclicks = 0;

	//
	selectin = frame;
	tmpimg = cvCloneImage(selectin);

	//
	cvShowImage(windowname, selectin);

	cvSetMouseCallback(windowname, mouse_callback, 0);

	//
	mouseeventprocessing = 1;
	stop = 0;

	while(!stop)
	{
		int key;

		//
		key = cvWaitKey(5);

		#define KEY_SPACE 32
		if(key == KEY_SPACE)
			stop = 1;
	}

	mouseeventprocessing = 0;

	if(nclicks == 2)
	{
		*pr = _r;
		*pc = _c;
		*ps = _s;
	}

	//
	cvReleaseImage(&tmpimg);
	tmpimg = 0;

	cvDestroyWindow(windowname);

	//
	return (nclicks==2);
}

/*
	
*/

int process_webcam_frames(IplImage** pimg, int* pr, int* pc, int* ps)
{
	CvCapture* capture;

	IplImage* frame;
	IplImage* framecopy;

	int stop, success;

	//
	cvNamedWindow(windowname, 1);

	// try to initialize video capture from the default webcam
	capture = cvCaptureFromCAM(0);
	///capture = cvCaptureFromAVI(path);

	if(!capture)
	{
		printf("cannot initialize video capture!\n");
		return 0;
	}

	// start the main loop in which we'll process webcam output
	framecopy = 0;

	stop = 0;
	success = 0;

	while(!stop)
	{
		int key;

		// wait ...
		key = cvWaitKey(5);

		// retrieve a pointer to the image acquired from the webcam
		if(!cvGrabFrame(capture))
			break;

		frame = cvRetrieveFrame(capture, 1);

		// we terminate the loop if we don't get any data from the webcam or the user has pressed 'q'
		if(!frame || key=='q')
			stop = 1;
		else
		{
			// we mustn't tamper with internal OpenCV buffers and that's the reason why we're making a copy of the current frame
			if(!framecopy)
				framecopy = cvCreateImage(cvSize(frame->width, frame->height), frame->depth, frame->nChannels);

			cvCopy(frame, framecopy, 0);

			// webcam outputs mirrored frames (at least on my machines); you can safely comment out this line if you find it unnecessary
			cvFlip(framecopy, framecopy, 1);

			//
			if(key == KEY_SPACE)
			{
				if( select_region(framecopy, pr, pc, ps) )
				{
					*pimg = cvCreateImage(cvSize(framecopy->width, framecopy->height), framecopy->depth, framecopy->nChannels);

					cvCopy(framecopy, *pimg, 0);

					success = 1;
					stop = 1;
				}
				else
					success = 0;
			}

			// display the image to the user
			cvShowImage(windowname, framecopy);
		}
	}

	// cleanup
	cvReleaseImage(&framecopy);
	cvReleaseCapture(&capture);
	cvDestroyWindow(windowname);

	//
	return success;
}

/*
	
*/

int main(int argc, char* argv[])
{
	int i;

	IplImage* gray;
	IplImage* rgb;

	int r, c, s, nwarps;

	char tag[1024], folder[1024];

	IplImage** warps;

	//
	if(argc == 8)
	{

		//
		rgb = cvLoadImage(argv[1], CV_LOAD_IMAGE_COLOR);
		gray = cvLoadImage(argv[1], CV_LOAD_IMAGE_GRAYSCALE);
	
		if(!gray || !rgb)
			return 0;
	
		sscanf(argv[2], "%d", &r);
		sscanf(argv[3], "%d", &c);
		sscanf(argv[4], "%d", &s);
		sscanf(argv[5], "%d", &nwarps);
	
		sscanf(argv[6], "%s", tag);
		sscanf(argv[7], "%s", folder);
	}
	else if(argc == 5)
	{
		//
		rgb = cvLoadImage(argv[1], CV_LOAD_IMAGE_COLOR);
		gray = cvLoadImage(argv[1], CV_LOAD_IMAGE_GRAYSCALE);

		if(!gray || !rgb)
			return 0;

		if(!select_region(rgb, &r, &c, &s))
			return 0;

		sscanf(argv[2], "%d", &nwarps);

		sscanf(argv[3], "%s", tag);
		sscanf(argv[4], "%s", folder);
	}
	else if(argc == 4)
	{
		//
		if(!process_webcam_frames(&rgb, &r, &c, &s))
			return 0;

		gray = cvCreateImage(cvSize(rgb->width, rgb->height), IPL_DEPTH_8U, 1);

		cvCvtColor(rgb, gray, CV_BGR2GRAY);

		//
		sscanf(argv[1], "%d", &nwarps);

		sscanf(argv[2], "%s", tag);
		sscanf(argv[3], "%s", folder);
	}
	else
	{
		printf("/path/to/image r c s nwarps tag /output/folder\n");

		return 0;
	}

	//
	warps = (IplImage**)malloc(nwarps*sizeof(IplImage*));

	generate_warps(&gray, &r, &c, &s, 1, nwarps, warps);

	//
	for(i=0; i<nwarps; ++i)
	{
		char path[1024];

		//
		sprintf(path, "%s/%s-%06d.rid", folder, tag, i);

		//
		saverid(path, warps[i]->imageData, warps[i]->height, warps[i]->width, warps[i]->widthStep);

		//
		printf("%s-%06d.rid\n", tag, i);
	}

	//
	return 0;
}