#ifndef __MRI_SCANNER_H
#define __MRI_SCANNER_H

//==========================================================================
// MRI_SCANNER.H
// MRI_Scanner class.
// Inherits from:  
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/mriscanner.h,v 1.2 2008-11-06 10:58:23 rotor Exp $
 * $Log: mriscanner.h,v $
 * Revision 1.2  2008-11-06 10:58:23  rotor
 *  * fixed includes for iostream and friends
 *  * updated for new release (1.0.2)
 *
 * Revision 1.1  2003/05/30 16:43:11  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:57:09  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:13:13  rkwan
 * Release 2.5
 *
 * Revision 1.4  1996/01/17  18:01:40  rkwan
 * Update for tx/rx coil modelling.
 *
 * Revision 1.3  1995/12/22  20:16:12  rkwan
 * Update for percent coil and RF map features.
 *
 * Revision 1.2  1995/12/11  15:23:50  rkwan
 * Doc update.
 *
 *========================================================================*/

#include <iostream>
#include "phantom.h"
#include "rf_coil.h"

using namespace std;

#include <minc/mriimage.h>
#include <minc/mrimatrix.h>
#include <minc/omincfile.h>
#include <signal/pulseseq.h>

//--------------------------------------------------------------------------
// MRI_Scanner class
//--------------------------------------------------------------------------

class MRI_Scanner {
   public:
      MRI_Scanner();

      virtual ~MRI_Scanner();

      // --- Register scanner components --- //

      int attach(Phantom *phantom);
      int attach(RF_Coil *rf_coil);
      int apply(Quick_Sequence *pseq);
      int apply(Custom_Sequence *pseq);

      inline int has_attached_phantom(void) const;
      inline int has_attached_rf_coil(void) const;
      inline int has_applied_pulse_sequence(void) const;

      inline Phantom *get_attached_phantom(void) const;
      inline Image_Type get_image_type(void) const;

      inline int save_raw_data(void) const;
      inline void save_raw_data(int on);

      // --- Volume information convenience functions --- //

      inline int    is_same_slice_size_as(const MRI_Matrix& mat) const;
      inline int    get_nrows(void) const;
      inline int    get_ncols(void) const;
      inline int    get_nslices(void) const;
      inline double get_voxel_step(int n) const;
      inline double get_voxel_start(int n) const;
      inline double get_voxel_offset(int n) const;

      inline int    get_matrix_size(int n) const;
      inline int    get_recon_size(int n) const;
   
      // --- Convenience functions --- //

      void display_info(ostream& stream) const ;
      int  initialize_output(O_MINC_File& output, const char *argstring) const;
      int  set_output_icv(O_MINC_File& output) const;

      inline double get_signal_gain(void) const; 
      inline void   set_signal_gain(double gain);

      // --- Simulated data generation --- //

      void get_simulated_image_slice(int slice, MRI_Image& image);
      void get_raw_data_slice(int slice, Complex_Slice& raw_slice);

      void initialize_chirp_resample(void);
      void reconstruct_raw_data_slice(const Complex_Slice& raw_slice, 
                                      Complex_Slice& output_slice);
      void reconstruct_raw_data_slice(Complex_Slice& raw_slice);

      void get_real_image(const Complex_Slice& complex_slice, 
                          MRI_Image& image_slice);
      void get_imag_image(const Complex_Slice& complex_slice,
                          MRI_Image& image_slice);
      void get_abs_image(const Complex_Slice& complex_slice,
                          MRI_Image& image_slice);
      void get_angle_image(const Complex_Slice& complex_slice,
                          MRI_Image& image_slice);

   private:

      // --- Internal data structures --- //

      Pulse_Sequence *_current_pseq;
      Phantom        *_phantom;
      RF_Coil        *_rf_coil;

      Volume_Info    _volume_info;
      double         _signal_gain;
      double         _voxel_offset[3];
    
      int            _save_raw_data;
 
      // --- Internal member functions --- // 

      void _update_volume_info(Phantom *phantom, Pulse_Sequence *pseq);

};

//--------------------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// MRI_Scanner::has_attached_phantom
// Returns TRUE if a phantom has been attached.
//--------------------------------------------------------------------------

inline 
int MRI_Scanner::has_attached_phantom(void) const {
   return (_phantom != NULL);
}

//--------------------------------------------------------------------------
// MRI_Scanner::has_attached_rf_coil
// Returns TRUE if an RF_Coil has been attached.
//--------------------------------------------------------------------------

inline 
int MRI_Scanner::has_attached_rf_coil(void) const {
   return (_rf_coil != NULL);
}

//--------------------------------------------------------------------------
// MRI_Scanner::has_applied_pulse_sequence
// Returns TRUE if a pulse sequence has been applied.
//--------------------------------------------------------------------------

inline
int MRI_Scanner::has_applied_pulse_sequence(void) const {
   return (_current_pseq != NULL);
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_attached_phantom
// Returns a pointer to the currently attached phantom.
//--------------------------------------------------------------------------

inline
Phantom *MRI_Scanner::get_attached_phantom(void) const {
   return _phantom;
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_image_type
// Returns the output type of the reconstructed image.
//--------------------------------------------------------------------------

inline
Image_Type MRI_Scanner::get_image_type(void) const {
   return _current_pseq->get_image_type();
}

//--------------------------------------------------------------------------
// MRI_Scanner::save_raw_data
// Returns TRUE if the save raw data option is turned on.
//--------------------------------------------------------------------------

inline
int MRI_Scanner::save_raw_data(void) const {
   return _save_raw_data;
}

//--------------------------------------------------------------------------
// MRI_Scanner::save_raw_data
// Turns on/off save raw data option.  on = TRUE, off = FALSE.
//--------------------------------------------------------------------------

inline
void MRI_Scanner::save_raw_data(int on) {
   _save_raw_data = on;
}

//--------------------------------------------------------------------------
// MRI_Scanner::is_same_slice_size_as
//--------------------------------------------------------------------------

inline 
int    MRI_Scanner::is_same_slice_size_as(const MRI_Matrix& mat) const {
   return ((mat.get_nrows() == _volume_info.length[ROW]) &&
           (mat.get_ncols() == _volume_info.length[COLUMN])); 
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_nrows
//--------------------------------------------------------------------------

inline 
int    MRI_Scanner::get_nrows(void) const {
   return _volume_info.length[ROW];
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_ncols 
//--------------------------------------------------------------------------

inline 
int    MRI_Scanner::get_ncols(void) const {
   return _volume_info.length[COLUMN];
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_nslices
//--------------------------------------------------------------------------

inline 
int    MRI_Scanner::get_nslices(void) const {
   return _volume_info.length[SLICE];
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_voxel_step
//--------------------------------------------------------------------------

inline 
double MRI_Scanner::get_voxel_step(int n) const {
   return _volume_info.step[n];
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_voxel_start
//--------------------------------------------------------------------------

inline 
double MRI_Scanner::get_voxel_start(int n) const {
   return _volume_info.start[n];
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_voxel_offset
//--------------------------------------------------------------------------

inline
double MRI_Scanner::get_voxel_offset(int n) const {
   return _voxel_offset[n];
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_matrix_size
//--------------------------------------------------------------------------

inline
int MRI_Scanner::get_matrix_size(int n) const {
   return _current_pseq->get_matrix_size(n);
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_recon_size
//--------------------------------------------------------------------------

inline
int MRI_Scanner::get_recon_size(int n) const {
   return _current_pseq->get_reconstruction_size(n);
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_signal_gain
// Returns the MRI_Scanner output image signal gain.
//--------------------------------------------------------------------------

inline 
double MRI_Scanner::get_signal_gain(void) const {
   return _signal_gain;
}

//--------------------------------------------------------------------------
// MRI_Scanner::set_signal_gain
// Sets the MRI_Scanner output image signal gain.
//--------------------------------------------------------------------------

inline 
void   MRI_Scanner::set_signal_gain(double gain) {
   _signal_gain = gain;
}

#endif

