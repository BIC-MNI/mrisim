//==========================================================================
// MRISCANNER.CXX
// MRI_Scanner class.
// Inherits from:  
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/mriscanner.cxx,v 1.2 2008-11-06 10:58:23 rotor Exp $
// $Log: mriscanner.cxx,v $
// Revision 1.2  2008-11-06 10:58:23  rotor
//  * fixed includes for iostream and friends
//  * updated for new release (1.0.2)
//
// Revision 1.1  2003/05/30 16:43:11  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  15:57:18  rkwan
// Release 3.1 update.
//
// Revision 2.6  1996/05/29  19:14:59  rkwan
// GNU warning fix.
//
// Revision 2.5  1996/05/29  16:13:22  rkwan
// Release 2.5
//
// Revision 1.5  1996/02/20  17:14:41  rkwan
// Fix volume start to compensate for position in the centre of voxels.
//
// Revision 1.4  1996/01/17  18:01:50  rkwan
// Update for tx/rx coil modelling.
//
// Revision 1.3  1995/12/22  20:15:57  rkwan
// Update for percent coil and RF map features.
//
// Revision 1.2  1995/12/11  15:24:02  rkwan
// Doc update.
//
//==========================================================================

#include <iostream>
#include <mrisim/mrisim.h>
#include "phantom.h"
#include "../minc/mriimage.h"
#include "../signal/pulseseq.h"

#include "mriscanner.h"

//--------------------------------------------------------------------------
// MRI_Scanner constructors
//--------------------------------------------------------------------------

MRI_Scanner::MRI_Scanner() {
   _current_pseq = (Pulse_Sequence *)NULL;
   _phantom      = (Phantom *)NULL;
   _rf_coil      = (RF_Coil *)NULL;

   _signal_gain  = 1.0;
}

//--------------------------------------------------------------------------
// MRI_Scanner destructor
//--------------------------------------------------------------------------

MRI_Scanner::~MRI_Scanner() {
   if (_phantom != NULL) delete _phantom;
   if (_rf_coil != NULL) delete _rf_coil;
}

//--------------------------------------------------------------------------
// MRI_Scanner::attach
// Attaches phantom and coil models to a scanner.
// Overloaded for the required model type.
//--------------------------------------------------------------------------

int MRI_Scanner::attach(Phantom *phantom) {

   int status;

   // Save a pointer to the scanner's phantom model
   if (phantom == NULL) {
      status = FALSE;
   } else {
      _phantom = phantom;
      status = TRUE;
   }
   return status;

}

int MRI_Scanner::attach(RF_Coil *rf_coil) {

   int status;

   // Save a pointer to the scanner's coil model
   // and initialize the coil model

   if (rf_coil == NULL) {
      status = FALSE;
   } else {
      _rf_coil = rf_coil;
      status = TRUE;
   }

   //_rf_coil->initialize_coil();
   return status;
}

//--------------------------------------------------------------------------
// MRI_Scanner::apply
// Applies the current pulse sequence to the phantom model.
//--------------------------------------------------------------------------

int MRI_Scanner::apply(Quick_Sequence *pseq) {

#ifdef DEBUG
   assert(pseq != NULL);
#endif
   
   _current_pseq = (Pulse_Sequence *)pseq;

   if (!has_attached_phantom()) {
      cerr << "MRI_Scanner:  No attached phantom." << endl;
      return FALSE;
   }

   _update_volume_info(_phantom, _current_pseq);
 
   _phantom->apply_pulse_sequence(pseq);
   _rf_coil->_compute_variance(*_phantom);

   return TRUE;
}

int MRI_Scanner::apply(Custom_Sequence *pseq) {

#ifdef DEBUG
   assert(pseq != NULL);
#endif

   _current_pseq = (Pulse_Sequence *)pseq;

   if (!has_attached_phantom()) {
      cerr << "MRI_Scanner:  No attached phantom." << endl;
      return FALSE;
   }

   _update_volume_info(_phantom, _current_pseq);

   _phantom->find_steady_state(pseq);
   _rf_coil->_compute_variance(*_phantom);

   return TRUE;
}

//--------------------------------------------------------------------------
// MRI_Scanner::display_info
// Writes information about the scanner models to an output stream.
//--------------------------------------------------------------------------

void MRI_Scanner::display_info(ostream& stream) const {

   stream << endl;

   stream << "MRI_Scanner Information:" << endl;
   stream << "------------------------" << endl << endl;
   stream << "signal gain: " << get_signal_gain() << endl << endl;

   if (this->has_applied_pulse_sequence()) {
      _current_pseq->display_info(stream);
   } else {
      stream << "MRI_Scanner WARNING:  No applied pulse sequence. " << endl;
   }

   if (this->has_attached_phantom()) {
      _phantom->display_volume_info(stream);
      _phantom->display_tissue_info(stream);
   } else {
      stream << "MRI_Scanner WARNING: No phantom attached." << endl;
   }

   if (this->has_attached_rf_coil()) {
      _rf_coil->display_coil_info(stream);
   } else {
      stream << "MRI_Scanner WARNING: No RF Coil attached." << endl;
   }

}

//--------------------------------------------------------------------------
// MRI_Scanner::initialize_output
// Sets up an output MINC file's volume information, creates and attaches
// an ICV to the output.
//--------------------------------------------------------------------------

int MRI_Scanner::initialize_output(O_MINC_File& output, 
                                   const char *argstring) const {

#ifdef DEBUG
   // Output should not have been previously initialized.
   assert(!output.is_good());
#endif

   if (!has_attached_phantom()) {
      cerr << "MRI_Scanner::initialize_output:  no phantom attached." 
           << endl;
      return FALSE;
   }

   _phantom->set_output_volume_info(output, _volume_info, argstring);
 
   if (!set_output_icv(output)) {
      return FALSE;
   }
   output.attach_icv();

   return TRUE;

}
   
//--------------------------------------------------------------------------
// MRI_Scanner::set_output_icv
// Initializes an output icv for the simulated images.
//--------------------------------------------------------------------------

int MRI_Scanner::set_output_icv(O_MINC_File& output) const {

   int status;

   if (!has_applied_pulse_sequence()) {
      cerr << "MRI_Scanner::set_output_icv: no applied pulse sequence."
           << endl;
      status = FALSE;
   } else {
      output.set_icv_property(MI_ICV_TYPE, (int)_volume_info.datatype);
      output.set_icv_property(MI_ICV_SIGN, (char *)_volume_info.signtype);
      output.set_icv_property(MI_ICV_VALID_MIN, _volume_info.valid_range[0]);
      output.set_icv_property(MI_ICV_VALID_MAX, _volume_info.valid_range[1]);
      output.set_icv_property(MI_ICV_DO_NORM, FALSE);

      status = TRUE;
   }

   return status;
}

//--------------------------------------------------------------------------
// MRI_Scanner::get_simulated_image_slice
// Computes a simulated image according to the scanner models.
//--------------------------------------------------------------------------

void MRI_Scanner::get_simulated_image_slice(int slice, MRI_Image& image_slice) {

   const double slice_thickness  = _current_pseq->get_slice_thickness();
   const double slice_separation = _current_pseq->get_voxel_step(SLICE);
   const double z_centre         = slice * slice_separation + 
                                   get_voxel_offset(SLICE);

   Real_Slice   phantom_slice(_phantom->get_nrows(), _phantom->get_ncols());
   Real_Slice   noisy_slice(image_slice.get_nrows(), image_slice.get_ncols());

   // Get partial volume image
   _phantom->ideal_nn_slice_select(z_centre, slice_thickness, phantom_slice); 
   _phantom->compute_partial_volume(phantom_slice, 
                             get_voxel_step(ROW),   get_voxel_step(COLUMN),
                             get_voxel_offset(ROW), get_voxel_offset(COLUMN),
                             noisy_slice);

   // Add noise and scale the image

   switch(get_image_type()){
      case REAL_IMAGE:
      case IMAG_IMAGE:
         _rf_coil->add_noise_to_real_image(noisy_slice);
         break;
      case MOD_IMAGE:
         _rf_coil->add_noise_to_modulus_image(noisy_slice);
      default:
         break;
   }

   image_slice.set_real_minimum(noisy_slice.minimum());
   image_slice.set_real_maximum(noisy_slice.maximum());

   double value;
   unsigned int m,n;
   for (m=0; m<image_slice.get_nrows(); m++){
      for (n=0; n<image_slice.get_ncols(); n++){
         value = noisy_slice(m,n);
         image_slice.set_value(m, n, value);
      }
   }

}

//--------------------------------------------------------------------------
// MRI_Scanner::get_raw_slice
// Computes raw fourier data according to the scanner models.
//--------------------------------------------------------------------------

void MRI_Scanner::get_raw_data_slice(int slice, Complex_Slice& raw_slice) {

   const double slice_thickness  = _current_pseq->get_slice_thickness();
   const double slice_separation = _current_pseq->get_voxel_step(SLICE);
   const double z_centre         = slice * slice_separation + 
                                   get_voxel_offset(SLICE);

   Real_Slice   phantom_slice(_phantom->get_nrows(), _phantom->get_ncols());

   _phantom->ideal_lin_slice_select(z_centre, slice_thickness, phantom_slice);
   _phantom->generate_raw_data_slice(phantom_slice, raw_slice);

   _rf_coil->add_noise_to_raw_slice(raw_slice);

}

//--------------------------------------------------------------------------
// MRI_Scanner::initialize_chirp_resample
// Initializes the Chirp DFT filters need for Chirp resampling.
//--------------------------------------------------------------------------

void MRI_Scanner::initialize_chirp_resample(void) {

   // The 2-D Chirp DFT is performed by a row-column decomposition of
   // 1-D Chirp DFTs.  A 1-D Chirp is performed on each row of the matrix
   // followed by a 1-D Chirp on each column of the matrix.
   // Note:   a row is a 1-D slice taken along the COLUMN direction, and
   // a col is a 1-D slice taken along the ROW direction.   There are
   // thus get_nrows() rows each of length get_ncols() elements and
   // intersample distance get_voxel_step(COLUMN).

   const unsigned int out_row_length = _current_pseq->get_matrix_size(COLUMN);
   const unsigned int out_col_length = _current_pseq->get_matrix_size(ROW);

   const double out_row_fov = out_row_length * 
                              _current_pseq->get_voxel_step(COLUMN);
   const double out_col_fov = out_col_length * 
                              _current_pseq->get_voxel_step(ROW);

   _phantom->initialize_chirp(out_row_length, out_col_length,
                              out_row_fov,    out_col_fov);

}

//--------------------------------------------------------------------------
// MRI_Scanner::reconstruct_raw_data_slice
// Reconstructs an MR image from the complex raw data slice.
//--------------------------------------------------------------------------

void MRI_Scanner::reconstruct_raw_data_slice(const Complex_Slice& raw_slice,
                                             Complex_Slice& output_slice) {

   int m, n;
   for (m=0; m<output_slice.get_nrows(); m++){
      for (n=0; n<output_slice.get_ncols(); n++){
         output_slice.real(m,n) = raw_slice.real(m,n);
         output_slice.imag(m,n) = raw_slice.imag(m,n);
      }
   }
   output_slice.fftshift();
   output_slice.iFFT2();

}

void MRI_Scanner::reconstruct_raw_data_slice(Complex_Slice &raw_slice) {

   raw_slice.fftshift();
   raw_slice.iFFT2();

}

//--------------------------------------------------------------------------
// MRI_Scanner::get_real_image
// Extracts a real image from a complex slice.
//--------------------------------------------------------------------------

void MRI_Scanner::get_real_image(const Complex_Slice& complex_slice,
                                 MRI_Image& image_slice) {

   double min, max;
   complex_slice.get_real_min_max(min, max);
   image_slice.set_real_minimum(min);
   image_slice.set_real_maximum(max);

   unsigned int m,n;
   for (m=0; m<image_slice.get_nrows(); m++){
      for (n=0; n<image_slice.get_ncols(); n++){
         image_slice.set_value(m, n, complex_slice.real(m,n));
      }
   }

}

//--------------------------------------------------------------------------
// MRI_Scanner::get_imag_image
// Extracts an imaginary image from a complex slice.
//--------------------------------------------------------------------------

void MRI_Scanner::get_imag_image(const Complex_Slice& complex_slice,
                                 MRI_Image& image_slice) {

   double min, max;
   complex_slice.get_imag_min_max(min, max);
   image_slice.set_real_minimum(min);
   image_slice.set_real_maximum(max);

   unsigned int m,n;
   for (m=0; m<image_slice.get_nrows(); m++){
      for (n=0; n<image_slice.get_ncols(); n++){
         image_slice.set_value(m, n, complex_slice.imag(m,n));
      }
   }

}

//--------------------------------------------------------------------------
// MRI_Scanner::get_abs_image
// Extracts a magnitude image from a complex slice.
//--------------------------------------------------------------------------
void MRI_Scanner::get_abs_image(const Complex_Slice& complex_slice,
                                MRI_Image& image_slice) {

   double min, max;
   complex_slice.get_abs_min_max(min, max);
   image_slice.set_real_minimum(min);
   image_slice.set_real_maximum(max);

   double value;
   unsigned int m,n;
   for (m=0; m<image_slice.get_nrows(); m++){
      for (n=0; n<image_slice.get_ncols(); n++){
         value = hypot(complex_slice.real(m,n), complex_slice.imag(m,n));
         image_slice.set_value(m, n, value);
      }
   }

}

//--------------------------------------------------------------------------
// MRI_Scanner::get_angle_image
// Extracts a phase image from a complex slice.
//--------------------------------------------------------------------------

void MRI_Scanner::get_angle_image(const Complex_Slice& complex_slice,
                                  MRI_Image& image_slice) {

   double min, max;
   complex_slice.get_angle_min_max(min, max);
   image_slice.set_real_minimum(min);
   image_slice.set_real_maximum(max);

   double value;
   unsigned int m,n;
   for (m=0; m<image_slice.get_nrows(); m++){
      for (n=0; n<image_slice.get_ncols(); n++){
         value = atan2(complex_slice.imag(m,n), complex_slice.real(m,n));
         image_slice.set_value(m, n, value);
      }
   }

}

//--------------------------------------------------------------------------
// MRI_Scanner::_update_volume_info
// Updates the scanner's volume information according to the Phantom and
// current Pulse_Sequence objects used.
//--------------------------------------------------------------------------

void MRI_Scanner::_update_volume_info(Phantom *phantom, Pulse_Sequence *pseq) {

#ifdef DEBUG
   assert(phantom != NULL);
   assert(pseq != NULL);
#endif

   phantom->get_volume_info(_volume_info);

   _volume_info.length[SLICE]  = pseq->get_nslices();
   _volume_info.length[ROW]    = pseq->get_nrows();
   _volume_info.length[COLUMN] = pseq->get_ncols();

   _volume_info.step[SLICE]    = pseq->get_voxel_step(SLICE);
   _volume_info.step[ROW]      = pseq->get_voxel_step(ROW);   
   _volume_info.step[COLUMN]   = pseq->get_voxel_step(COLUMN);

   const double row_offset = (phantom->get_nrows() *
                              phantom->get_voxel_step(ROW) - 
                              pseq->get_fov(ROW))/2 +
                              pseq->get_voxel_offset(ROW);

   const double col_offset = (phantom->get_ncols() *
                              phantom->get_voxel_step(COLUMN) - 
                              pseq->get_fov(COLUMN))/2 +
                              pseq->get_voxel_offset(COLUMN);

   _voxel_offset[SLICE]  = pseq->get_voxel_offset(SLICE);
   _voxel_offset[ROW]    = row_offset;
   _voxel_offset[COLUMN] = col_offset;

   _volume_info.start[SLICE]   += pseq->get_voxel_offset(SLICE);
   _volume_info.start[ROW]     += row_offset;
   _volume_info.start[COLUMN]  += col_offset;

   _volume_info.datatype       = NC_SHORT;

   strcpy(_volume_info.signtype, MI_SIGNED);

   _volume_info.valid_range[0] = 0.0;
   _volume_info.valid_range[1] = 4095.0;

}
