/********************************************************************************
 * FOURN.C
 * C support functions for FFTs
 * 
 * R. Kwan
 * August 31, 1995
 *
 *******************************************************************************/

#include <math.h>
#include <string.h>
#include "fourn.h"

#ifdef TRACE
#include <stdio.h>
#endif

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/*
 * fourn
 * n-dimensional FFT code.  Taken from Numerical Recipes.
 */
void fourn(double data[],unsigned long nn[],int ndim,int isign)
{
        int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	double tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

        ntot=1;
        for (idim=1;idim<=ndim;idim++)
           ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
                ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=wr*data[k2]-wi*data[k2+1];
						tempi=wr*data[k2+1]+wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}

/*
 * four1
 * 1D-FFT.
 * Taken from Numerical Recipes.
 */

void four1(double data[], int nn, int isign){
   int n,mmax,m,j,istep,i;
   double wtemp,wr,wpr,wpi,wi,theta;
   double tempr,tempi;

   n=nn << 1;
   j=1;
   for (i=1;i<n;i+=2) {
      if (j > i) {
         SWAP(data[j],data[i]);
         SWAP(data[j+1],data[i+1]);
      }
      m=n >> 1;
      while (m >= 2 && j > m) {
         j -= m;
         m >>= 1;
      }
      j += m;
   }
   mmax=2;
   while (n > mmax) {
      istep=2*mmax;
      theta=6.28318530717959/(isign*mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=1;m<mmax;m+=2) {
         for (i=m;i<=n;i+=istep) {
            j=i+mmax;
            tempr=wr*data[j]-wi*data[j+1];
            tempi=wr*data[j+1]+wi*data[j];
            data[j]=data[i]-tempr;
            data[j+1]=data[i+1]-tempi;
            data[i] += tempr;
            data[i+1] += tempi;
         }
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
      } 
      mmax=istep;
   }
}
/*
 * fournf
 * n-dimensional FFT code.  Taken from Numerical Recipes.
 */
void fournf(float data[],unsigned long nn[],int ndim,int isign)
{
        int idim;
	unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
	unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
	float  tempi,tempr;
	double theta,wi,wpi,wpr,wr,wtemp;

        ntot=1;
        for (idim=1;idim<=ndim;idim++)
           ntot *= nn[idim];
	nprev=1;
	for (idim=ndim;idim>=1;idim--) {
		n=nn[idim];
		nrem=ntot/(n*nprev);
                ip1=nprev << 1;
		ip2=ip1*n;
		ip3=ip2*nrem;
		i2rev=1;
		for (i2=1;i2<=ip2;i2+=ip1) {
			if (i2 < i2rev) {
				for (i1=i2;i1<=i2+ip1-2;i1+=2) {
					for (i3=i1;i3<=ip3;i3+=ip2) {
						i3rev=i2rev+i3-i2;
						SWAP(data[i3],data[i3rev]);
						SWAP(data[i3+1],data[i3rev+1]);
					}
				}
			}
			ibit=ip2 >> 1;
			while (ibit >= ip1 && i2rev > ibit) {
				i2rev -= ibit;
				ibit >>= 1;
			}
			i2rev += ibit;
		}
		ifp1=ip1;
		while (ifp1 < ip2) {
			ifp2=ifp1 << 1;
			theta=isign*6.28318530717959/(ifp2/ip1);
			wtemp=sin(0.5*theta);
			wpr = -2.0*wtemp*wtemp;
			wpi=sin(theta);
			wr=1.0;
			wi=0.0;
			for (i3=1;i3<=ifp1;i3+=ip1) {
				for (i1=i3;i1<=i3+ip1-2;i1+=2) {
					for (i2=i1;i2<=ip3;i2+=ifp2) {
						k1=i2;
						k2=k1+ifp1;
						tempr=(float)wr*data[k2]-
                                                      (float)wi*data[k2+1];
						tempi=(float)wr*data[k2+1]+
                                                      (float)wi*data[k2];
						data[k2]=data[k1]-tempr;
						data[k2+1]=data[k1+1]-tempi;
						data[k1] += tempr;
						data[k1+1] += tempi;
					}
				}
				wr=(wtemp=wr)*wpr-wi*wpi+wr;
				wi=wi*wpr+wtemp*wpi+wi;
			}
			ifp1=ifp2;
		}
		nprev *= n;
	}
}

/*
 * four1f
 * 1D-FFT.
 * Taken from Numerical Recipes.
 */

void four1f(float data[], int nn, int isign){
   int n,mmax,m,j,istep,i;
   double wtemp,wr,wpr,wpi,wi,theta;
   float tempr,tempi;

   n=nn << 1;
   j=1;
   for (i=1;i<n;i+=2) {
      if (j > i) {
         SWAP(data[j],data[i]);
         SWAP(data[j+1],data[i+1]);
      }
      m=n >> 1;
      while (m >= 2 && j > m) {
         j -= m;
         m >>= 1;
      }
      j += m;
   }
   mmax=2;
   while (n > mmax) {
      istep=2*mmax;
      theta=6.28318530717959/(isign*mmax);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (m=1;m<mmax;m+=2) {
         for (i=m;i<=n;i+=istep) {
            j=i+mmax;
            tempr=(float)wr*data[j]-(float)wi*data[j+1];
            tempi=(float)wr*data[j+1]+(float)wi*data[j];
            data[j]=data[i]-tempr;
            data[j+1]=data[i+1]-tempi;
            data[i] += tempr;
            data[i+1] += tempi;
         }
         wr=(wtemp=wr)*wpr-wi*wpi+wr;
         wi=wi*wpr+wtemp*wpi+wi;
      } 
      mmax=istep;
   }
}

/*
 * dftgl
 * 1D-Goertzel DFT algorithm.
 * Computes DFT components starting at kstart and ending at kend.
 * Adapted from FORTRAN code in Proakis & Manolakis,
 * "Digital Signal Processing".
 */

void dftgl(double data[], int N, double dataf[], int L, int points[]){

   double wn = 2*M_PI/(double)N;
   double vr1, vr2, vi1, vi2;
   double phi, c, s, temp;
   int i,k;
  
   for(k=0; k<2*L; k+=2){
      vr1 = 0.0;
      vr2 = 0.0;
      vi1 = 0.0;
      vi2 = 0.0;
      phi = wn*(points[k/2]);
      c   = cos(phi);
      s   = sin(phi);
      for(i=0; i<2*N; i+=2){
         temp = vr1;
         vr1 = 2*c*vr1 - vr2 + data[i];
         vr2 = temp;
         temp = vi1;
         vi1 = 2*c*vi1 - vi2 + data[i+1];
         vi2 = temp;
      }
      dataf[k]   = c*vr1 - vr2 - s*vi1;
      dataf[k+1] = c*vi1 - vi2 + s*vr1;
   }
}

/*
 * chirpdft
 * Chirp DFT
 * Returns L points of the DFT along exp(jtheta)*exp(jphi)^k,
 * in the array y.
 * Reference:  Oppenheim, A.V. & R.W. Schafer, "Discrete-Time Signal
 *             Processing", Prentice-Hall.
 *             Proakis & Manolakis, "Digital Signal Processing",
 *             MacMillan.
 */

void chirpdft(double in[], int in_length, double out[], int out_length, 
              double w_initial, double w_step){

   /* Define filter lengths */
   int filt_length  = in_length+out_length-1;
   int chirp_length = ((in_length > out_length) ? 2*in_length-1 : filt_length);
   int fft_length   = _next_power_of_two(filt_length);
   int offset       = 2*(in_length-1);

   /* Allocate filters */
   double *chirp     = (double *)calloc(2*chirp_length, sizeof(double));
   double *chirp_fft = (double *)calloc(2*fft_length, sizeof(double));
   double *tmp       = (double *)calloc(2*fft_length, sizeof(double));

   /* Offset filters */
   double *tmp_no_alias     = &(tmp[offset]);
   double *chirp_zero_delay = &(chirp[offset]);

   /* Define temporaries */
   int n;
   double tempr, tempi, filtr, filti;

   if ((tmp==NULL) || (chirp==NULL) || (chirp_fft==NULL)){
      free(tmp); free(chirp); free(chirp_fft);
      return;
   }

   /* Compute chirp filter */
   for(n=0; n<2*chirp_length; n+=2){
      chirp[n]   =  cos(w_step*(n-offset)*(n-offset)/8.0);
      chirp[n+1] = -sin(w_step*(n-offset)*(n-offset)/8.0);
   }
   
   /* Compute FFT of chirp filter for convolution */
   for(n=0; n<2*fft_length; n+=2){
      chirp_fft[n]   =  chirp[n];
      chirp_fft[n+1] = -chirp[n+1];
   }
   four1(chirp_fft-1, fft_length, -1);

   /* Premultiply data */
   for(n=0; n<2*in_length; n+=2){
      tempr =  cos(w_initial*n/2.0);
      tempi = -sin(w_initial*n/2.0);
      filtr =  tempr*chirp_zero_delay[n] - tempi*chirp_zero_delay[n+1];
      filti =  tempi*chirp_zero_delay[n] + tempr*chirp_zero_delay[n+1];
      
      tmp[n]   =  filtr*in[n] - filti*in[n+1];
      tmp[n+1] =  filti*in[n] + filtr*in[n+1];
   }
   /* Zero pad to FFT length */
   for(n=2*in_length; n<2*fft_length; n++){
      tmp[n] = 0.0;
   }

   /* Perform fast convolution via FFT */
   four1(tmp-1, fft_length, -1);
   for(n=0; n<2*fft_length; n+=2){
      tempr    = tmp[n];
      tempi    = tmp[n+1];
      tmp[n]   = tempr*chirp_fft[n] - tempi*chirp_fft[n+1];
      tmp[n+1] = tempi*chirp_fft[n] + tempr*chirp_fft[n+1];

   }
   four1(tmp-1, fft_length, 1);

   /* Postmultiply data */
   for(n=0; n<2*out_length; n+=2){
      out[n]   =  ( tmp_no_alias[n]*chirp_zero_delay[n] - 
                    tmp_no_alias[n+1]*chirp_zero_delay[n+1])/fft_length;
      out[n+1] =  ( tmp_no_alias[n]*chirp_zero_delay[n+1] + 
                    tmp_no_alias[n+1]*chirp_zero_delay[n])/fft_length;
   }
   free(chirp); free(chirp_fft); free(tmp);
} 
 
/*
 * memswap
 * Swaps n bytes of memory from location s1 and s2.
 * Returns a pointer to s1 or NULL if an error has occurred.
 */

void *memswap(void *s1, void *s2, size_t n){

   char *tmp;

   if ((tmp = (char *)malloc(n)) == NULL){
      return NULL;
   }   
   (void)memcpy(tmp,s1,n);
   (void)memcpy(s1,s2,n);
   (void)memcpy(s2,tmp,n);

   free(tmp);
   return s1;

}

/*
 * _power_of_two
 * Returns 1 if n is a power of two and 0 otherwise.
 */

int _power_of_two(unsigned int n){
   int p = 1, result = 0;

   while ((p < n) && (result == 0)){
      if (p == n) {
         result = 1;
      } else {
         p *= 2;
      }
   }
   return result;
}

/*
 * _next_power_of_two
 * Rounds up to the next highest power of two
 */

int _next_power_of_two(unsigned int n){
   unsigned p = 1;
   while (p < n){
      p *= 2;
   }
   return p;
}

/*
 * fftshift
 * Swaps first and third quadrants and second and fourth quadrants, to move the
 * zeroth lag to the center of the spectrum.
 * Needs the number of rows and cols in the matrix as well as the size in
 * bytes of each matrix elements.   The matrix should be stored by row.
 */

void _fftshift(void *mat, unsigned int nrows, unsigned int ncols, size_t el_size){
   int   row, col;
   int   M, N;
   int   offset1, offset2;
   unsigned char *ptr, tmp;

   /* Compute half-way quadrant divisions
    * and array offsets between quadrants
    */
   M = (int)ceil(nrows/2);
   N = (int)ceil(ncols/2);

   offset1 = el_size*(M*ncols+N);  /* offset between quadrants 1 & 4 */
   offset2 = el_size*(M*ncols-N);  /* offset between quadrants 2 & 3 */

   for(ptr=(unsigned char *)mat, row=0; row<M; row++){
      /* Swap first and fourth quadrants */
      for(col=0; col<el_size*N; col++, ptr++){
         tmp = *ptr;
         *ptr = *(ptr+offset1);
         *(ptr+offset1) = tmp;
      }
      /* Swap second and third quadrants */
      for(col=el_size*N; col<el_size*ncols; col++, ptr++){
         tmp = *ptr;
         *ptr = *(ptr+offset2);
         *(ptr+offset2) = tmp;
      }
   }
}  

void _fftshift_1d(void *mat, unsigned int N, size_t el_size){
   int offset = el_size*(N/2);
   memswap(mat, (char *)mat+offset, offset);
}

#undef SWAP
