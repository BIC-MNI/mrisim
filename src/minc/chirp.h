#ifndef __CHIRP_H
#define __CHIRP_H

//===========================================================================
// CHIRP.H
// Chirp Fourier Transform Algorithm.
// Inherits from:
// Base class to:
// R. Kwan
// June 20, 1996
// 
// (C) Copyright 1996 by R.Kwan
//===========================================================================

/*===========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/chirp.h,v 1.1 2003-05-30 16:43:08 bert Exp $
 * $Log: chirp.h,v $
 * Revision 1.1  2003-05-30 16:43:08  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:45:09  rkwan
 * Release 3.1 - initial release of chirp DFT resampling.
 *
 *=========================================================================*/

#include "fourn.h"

//---------------------------------------------------------------------------
// Chirp_Algorithm
// Encapsulates the Chirp Fourier Transform Algorithm.
//---------------------------------------------------------------------------

class Chirp_Algorithm {
   public:
      Chirp_Algorithm(int in_length, int out_length,
                      double w_initial, double w_step);
      virtual ~Chirp_Algorithm();

      // --- Access functions --- //

      inline int get_input_length(void) const;
      inline int get_output_length(void) const;
      inline double get_initial_freq(void) const;
      inline double get_step_freq(void) const;

      // --- Apply the chirp algorithm --- //

      void apply(int complex_input,
                 const float in[],  unsigned int in_stride,
                 float out[],       unsigned int out_stride);

   private:
      int _in_length;                   // input vector length
      int _out_length;                  // output vector length
      int _fft_length;                  // length for fast convolution

      double _w_initial;                // initial sample frequency
      double _w_step;                   // frequency step between samples

      double       *_chirp;             // the Chirp filter
      double       *_chirp_fft;         // FFT of the Chirp filter for conv.
      double       *_tmp;               // work space
      double       *_prefilter;         // the Chirp data pre-filter
      double       *_chirp_delay_zero;  // pointer to zero delay chirp filter
      double       *_tmp_no_alias;      // pointer to no alias section of conv.

};

//---------------------------------------------------------------------------
// Inline member functions
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Chirp_Algorithm::get_input_length
// Returns the input vector length of the Chirp DFT.
//---------------------------------------------------------------------------

inline 
int Chirp_Algorithm::get_input_length(void) const {
   return _in_length;
}

//---------------------------------------------------------------------------
// Chirp_Algorithm::get_output_length
// Returns the output vector length of the Chirp DFT.
//---------------------------------------------------------------------------

inline 
int Chirp_Algorithm::get_output_length(void) const {
   return _out_length;
}

//---------------------------------------------------------------------------
// Chirp_Algorithm::get_freq_initial
// Returns the frequency of the initial DFT sample.
//---------------------------------------------------------------------------

inline 
double Chirp_Algorithm::get_initial_freq(void) const {
   return _w_initial;
}

//---------------------------------------------------------------------------
// Chirp_Algorithm::get_freq_step
// Returns the frequency step of DFT samples.
//---------------------------------------------------------------------------

inline 
double Chirp_Algorithm::get_step_freq(void) const {
   return _w_step;
}

#endif
