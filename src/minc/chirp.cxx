//===========================================================================
// CHIRP.CXX
// 
// R.Kwan
// June 20, 1996
// 
// (C) Copyright 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/chirp.cxx,v 1.1 2003-05-30 16:43:08 bert Exp $
// $Log: chirp.cxx,v $
// Revision 1.1  2003-05-30 16:43:08  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  15:46:53  rkwan
// Release 3.1 - initial release of chirp DFT resampling.
//
//===========================================================================

#include <math.h>
#include <string.h>
#include "chirp.h"

//===========================================================================
// Chirp_Algorithm
//
// Encapsulates the Chirp Fourier Algorithm.
// Accepts an input vector of length in_length, and returns out_length
// points (k = 0, 1, ... , out_length-1) of the DFT along 
// exp(-j*w_initial)*exp(-j*k*w_step).
//
// Reference:  Oppenheim, A.V. & R.W. Schafer, "Discrete-Time Signal
//             Processing", Prentice-Hall.
//             Proakis & Manolakis, "Digital Signal Processing",
//             MacMillan.
//===========================================================================

//---------------------------------------------------------------------------
// Chirp_Algorithm constructor
//---------------------------------------------------------------------------

Chirp_Algorithm::Chirp_Algorithm(int in_length, 
                                 int out_length,
                                 double w_initial,
                                 double w_step) 
   : _in_length(in_length), _out_length(out_length), 
     _w_initial(w_initial), _w_step(w_step) {

#ifdef DEBUG
   assert(w_initial >= -M_PI);
   assert(w_initial <= M_PI);
   assert(w_initial + out_length*w_step >= w_initial);
   assert(w_initial + out_length*w_step <= M_PI);
#endif

   // --- Compute filter lengths --- //

   int chirp_length = ((in_length > out_length) ?
                        2*in_length-1 : in_length+out_length-1);
   int filt_length  = (in_length+out_length-1);

   _fft_length = _next_power_of_two(filt_length);

   const int offset = 2*(_in_length-1);

   // --- Allocate filters --- //
   // Filters are complex arrays arranged as 
   // {real0, imag0, real1, imag1, ... } 

   _chirp     = new double[2*chirp_length];
   _chirp_fft = new double[2*_fft_length];
   _tmp       = new double[2*_fft_length];
   _prefilter = new double[2*_in_length];

   _chirp_delay_zero = &(_chirp[offset]);
   _tmp_no_alias     = &(_tmp[offset]);

   // --- Compute chirp filter --- //

   // h[n] = exp(j*w_step)^(n^2/2) 
   // Chirp filter = h[n] n = -(in_length-1), ..., out_length-1 

   int n;
   double angle;

   for(n=0; n<2*chirp_length; n+=2){
      angle       =  w_step*(double)((n-offset)*(n-offset))/8.0;
      _chirp[n]   =  cos(angle);
      _chirp[n+1] = -sin(angle);
   }

   // --- Compute chirp filter FFT --- //

   for(n=0; n<2*filt_length; n+=2){
      _chirp_fft[n]   =  _chirp[n];
      _chirp_fft[n+1] = -_chirp[n+1];
   }
   for(n=2*filt_length; n<2*_fft_length; n++){
      _chirp_fft[n] = 0.0;
   }
   four1( _chirp_fft-1, _fft_length, -1);

   // --- Compute pre-filter --- //

   // Prefilter p[n] = exp(-j*w_initial*n) * h[n]^(-1)
   // n = 0, 1, ..., _in_length-1

   double tempr, tempi, filtr, filti;

   for(n=0; n<2*_in_length; n+=2){
      angle =  w_initial*(double)n/2.0;
      tempr =  cos(angle);
      tempi = -sin(angle); 
      filtr =  _chirp_delay_zero[n];
      filti =  _chirp_delay_zero[n+1];

      _prefilter[n]   = tempr * filtr - tempi * filti;
      _prefilter[n+1] = tempr * filti + tempi * filtr;
   }

}

//---------------------------------------------------------------------------
// Chirp_Algorithm destructor
//---------------------------------------------------------------------------

Chirp_Algorithm::~Chirp_Algorithm() {
   delete[] _chirp;
   delete[] _chirp_fft;
   delete[] _tmp;
   delete[] _prefilter;
}

//---------------------------------------------------------------------------
// Chirp_Algorithm::apply
//---------------------------------------------------------------------------

void Chirp_Algorithm::apply(int complex_input,
                            const float in[],  unsigned int in_stride,
                            float out[],       unsigned int out_stride) {

   int n, m;
   double tempr, tempi, filtr, filti;

   // --- Premultiply by chirp filter --- //

   if (complex_input) {

      for(n=0, m=0; n<2*_in_length; n+=2, m+=(2*in_stride)){
         tempr = in[m]; 
         tempi = in[m+1];
         filtr = _prefilter[n];
         filti = _prefilter[n+1];
 
         _tmp[n]   = tempr * filtr - tempi * filti;
         _tmp[n+1] = tempr * filti + tempi * filtr;
      }

   } else {  // real input

      for(n=0, m=0; n<2*_in_length; n+=2, m+=in_stride){
         tempr = in[m];

         _tmp[n]   = _prefilter[n]   * tempr;
         _tmp[n+1] = _prefilter[n+1] * tempr;
      }
   
   }

   // --- Zero-pad the rest of the pre-filtered input vector --- //

   for(n=2*_in_length; n<2*_fft_length; n++){
      _tmp[n] = 0.0;   
   }

   // --- Use FFT for fast convolution --- //

   four1(_tmp-1, _fft_length, -1);
   for(n=0; n<2*_fft_length; n+=2){
      tempr = _tmp[n];  
      tempi = _tmp[n+1];
      filtr = _chirp_fft[n];
      filti = _chirp_fft[n+1];

      _tmp[n]   = tempr * filtr - tempi * filti;
      _tmp[n+1] = tempr * filti + tempi * filtr;
   }
   four1(_tmp-1, _fft_length, 1);

   // --- Postmultiply data --- //
   // Note: postmultiplication is by 1./_chirp, which in
   // the Fourier case simplifies to the complex conjugate.

   for(n=0, m=0; n<2*_out_length; n+=2, m+=(2*out_stride)){
      tempr = _tmp_no_alias[n];
      tempi = _tmp_no_alias[n+1];
      filtr = _chirp_delay_zero[n];
      filti = _chirp_delay_zero[n+1];
     
      out[m]   = (tempr * filtr - tempi * filti) / _fft_length;
      out[m+1] = (tempr * filti + tempi * filtr) / _fft_length;
   }
                             
}
