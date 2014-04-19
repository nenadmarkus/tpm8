#ifndef ABS
#define ABS(x) ((x)>0?(x):(-(x)))
#endif

#define NORMALIZATION 255

int bintest(int32_t tcode, int threshold, int r, int c, int s, uint8_t pixels[], int nrows, int ncols, int ldim)
{
	int r1, c1, r2, c2;

	int8_t* ptr;

	//
	ptr = (int8_t*)&tcode;

	//
	r1 = (NORMALIZATION*r + ptr[0]*s)/NORMALIZATION;
	c1 = (NORMALIZATION*c + ptr[1]*s)/NORMALIZATION;

	r2 = (NORMALIZATION*r + ptr[2]*s)/NORMALIZATION;
	c2 = (NORMALIZATION*c + ptr[3]*s)/NORMALIZATION;

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