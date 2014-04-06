/*
	
*/

#define MAX_DIVS 16
#define N_BINS 256

void RCLAHEM(unsigned char imap[], unsigned char img[], int i0, int j0, int i1, int j1, int ldim, int clim)
{	
	//
	double p[N_BINS];
	double P[N_BINS];

	int i, j, k;
	int I, J;
	
	//
	I = i1 - i0 + 1;
	J = j1 - j0 + 1;

	//
	for(i=0; i<N_BINS; ++i)
		p[i] = 0.0;

	// compute histogram
	for(i=i0; i<=i1; ++i)
		for(j=j0; j<=j1; ++j)
		{
			k = img[i*ldim + j];

			p[k] = p[k] + 1.0/(I*J);
		}
		
	// clip the histogram (ideally, we should do a few iterations of this)
	for(k=0; k<N_BINS; ++k)
	{
		if(p[k] >= (double)clim/N_BINS)
		{
			//
			double d;
		
			//
			d = p[k] - (double)clim/N_BINS;
		
			//
			p[k] = (double)clim/N_BINS;
			
			// redistribute d
			for(i=0; i<N_BINS; ++i)
				p[i] += d/N_BINS;
		}
	}

	// compute cumulative histogram
	P[0] = p[0];
	for(i=1; i<N_BINS; ++i)
		P[i] = P[i-1] + p[i];
	
	// compute intensity map
	for(k=0; k<N_BINS; ++k)
	{
		imap[k] = (N_BINS-1)*P[k];
	}
}

void CLAHE(unsigned char out[], unsigned char in[], int nrows, int ncols, int ldim, int divr, int divc, int clim)
{
		//
	unsigned char imaps[MAX_DIVS][MAX_DIVS][N_BINS];
	
	int ics[MAX_DIVS], jcs[MAX_DIVS];
	
	int i, j, k, l, i0, j0, i1, j1, I, J;
	
	unsigned char v00, v01, v10, v11;
	
	//
	i0 = 0;
	j0 = 0;
	
	I = nrows;
	J = ncols;
	
	for(i=0; i<divr; ++i)
	{
		for(j=0; j<divc; ++j)
		{	
			//
			i0 = i*I/divr;
			j0 = j*J/divc;
			
			i1 = (i+1)*I/divr;
			j1 = (j+1)*J/divc;
			
			if(i1>=I)
				i1 = I - 1;
				
			if(j1 >= J)
				j1 = J - 1;
			
			//
			RCLAHEM(imaps[i][j], in, i0, j0, i1, j1, ldim, clim);
			
			//
			ics[i] = (i0 + i1)/2;
			jcs[j] = (j0 + j1)/2;
		}
	}
		
	// SPECIAL CASE: image corners
	for(i=0; i<ics[0]; ++i)
	{
		for(j=0; j<jcs[0]; ++j)
			out[i*ldim + j] = imaps[0][0][ in[i*ldim + j] ];
	
		for(j=jcs[divc-1]; j<J; ++j)
			out[i*ldim + j] = imaps[0][divc-1][ in[i*ldim + j] ];
	}
	
	for(i=ics[divr-1]; i<I; ++i)
	{
		for(j=0; j<jcs[0]; ++j)
			out[i*ldim + j] = imaps[divr-1][0][ in[i*ldim + j] ];
		
		for(j=jcs[divc-1]; j<J; ++j)	
			out[i*ldim + j] = imaps[divr-1][divc-1][ in[i*ldim + j] ];
	}
		
	// SPECIAL CASE: image boundaries
	for(k=0; k<divr-1; ++k)
	{
		for(i=ics[k]; i<ics[k+1]; ++i)
		{
			for(j=0; j<jcs[0]; ++j)
			{
				v00 = imaps[k+0][0][ in[i*ldim + j] ];
				v10 = imaps[k+1][0][ in[i*ldim + j] ];
				
				out[i*ldim + j] = ( (ics[k+1]-i)*v00 + (i-ics[k])*v10 )/(ics[k+1] - ics[k]);
			}
			
			for(j=jcs[divc-1]; j<J; ++j)
			{
				v01 = imaps[k+0][divc-1][ in[i*ldim + j] ];
				v11 = imaps[k+1][divc-1][ in[i*ldim + j] ];
				
				out[i*ldim + j] = ( (ics[k+1]-i)*v01 + (i-ics[k])*v11 )/(ics[k+1] - ics[k]);
			}
		}
	}
	
	for(k=0; k<divc-1; ++k)
		for(j=jcs[k]; j<jcs[k+1]; ++j)
		{
			for(i=0; i<ics[0]; ++i)
			{
				v00 = imaps[0][k+0][ in[i*ldim + j] ];
				v01 = imaps[0][k+1][ in[i*ldim + j] ];
				
				out[i*ldim + j] = ( (jcs[k+1]-j)*v00 + (j-jcs[k])*v01 )/(jcs[k+1] - jcs[k]);
			}
			
			for(i=ics[divr-1]; i<I; ++i)
			{
				v10 = imaps[divr-1][k+0][ in[i*ldim + j] ];
				v11 = imaps[divr-1][k+1][ in[i*ldim + j] ];
				
				out[i*ldim + j] = ( (jcs[k+1]-j)*v10 + (j-jcs[k])*v11 )/(jcs[k+1] - jcs[k]);
			}
		}
	
	//
	for(k=0; k<divr-1; ++k)
		for(l=0; l<divc-1; ++l)
			for(j=jcs[l]; j<jcs[l+1]; ++j)
				for(i=ics[k]; i<ics[k+1]; ++i)
				{
					unsigned char p;
					
					p = in[i*ldim + j];
				
					v00 = imaps[k+0][l+0][p];
					v01 = imaps[k+0][l+1][p];
					v10 = imaps[k+1][l+0][p];
					v11 = imaps[k+1][l+1][p];
					
					out[i*ldim + j] =
						(
							(ics[k+1]-i)*(jcs[l+1]-j)*v00 + (ics[k+1]-i)*(j-jcs[l])*v01 + (i-ics[k])*(jcs[l+1]-j)*v10 + (i-ics[k])*(j-jcs[l])*v11
						)/((ics[k+1] - ics[k])*(jcs[l+1] - jcs[l]));
				}
}