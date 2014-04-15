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