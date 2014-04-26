#ifndef ABS
#define ABS(x) ((x)>0?(x):(-(x)))
#endif

#define _FIXED_POINT_SCALE_ 256
#define _SQR_FIXED_POINT_SCALE_ (_FIXED_POINT_SCALE_*_FIXED_POINT_SCALE_)

int bintest(int32_t tcode, int threshold, int* T, uint8_t* pixels, int nrows, int ncols, int ldim)
{
	int r1, c1, r2, c2;

	int8_t* p;

	//
	p = (int8_t*)&tcode;

	//
	r1 = (T[0]*p[0] + T[1]*p[1] + T[2])/_FIXED_POINT_SCALE_;
	c1 = (T[3]*p[0] + T[4]*p[1] + T[5])/_FIXED_POINT_SCALE_;

	r2 = (T[0]*p[2] + T[1]*p[3] + T[2])/_FIXED_POINT_SCALE_;
	c2 = (T[3]*p[2] + T[4]*p[3] + T[5])/_FIXED_POINT_SCALE_;

	//
#ifndef USE_RGB
	return ( ABS(pixels[r1*ldim+c1]-pixels[r2*ldim+c2]) > threshold );
#else
	return (ABS(pixels[r1*ldim+3*c1+0]-pixels[r2*ldim+3*c2+0])>threshold)||(ABS(pixels[r1*ldim+3*c1+1]-pixels[r2*ldim+3*c2+1])>threshold)||(ABS(pixels[r1*ldim+3*c1+2]-pixels[r2*ldim+3*c2+2])>threshold);
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