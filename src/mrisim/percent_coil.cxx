//==========================================================================
// PERCENT_COIL.CXX
// Percent SNR Coil class.
// Inherits from:  RF_Coil
// Base class to:
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/percent_coil.cxx,v 1.1 2003-05-30 16:43:11 bert Exp $
// $Log: percent_coil.cxx,v $
// Revision 1.1  2003-05-30 16:43:11  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  16:35:17  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:17:38  rkwan
// Release 2.5
//
// Revision 1.3  1995/12/22  20:15:20  rkwan
// Update for percent coil and RF map features.
//
// Revision 1.2  1995/12/11  15:29:20  rkwan
// Doc update.
//
//==========================================================================

#include "percent_coil.h"

#ifdef DEBUG
#include <assert.h>
#endif

//--------------------------------------------------------------------------
// Percent_Coil Class
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Percent_Coil Constructors
//--------------------------------------------------------------------------

Percent_Coil::Percent_Coil(long seed, const char *rx_map_file,
                           const char *tx_map_file) 
   : RF_Coil(seed, rx_map_file, tx_map_file) {
   _noise_percentage = 0;
}

Percent_Coil::Percent_Coil(double percent, 
                           double reference_thickness,
                           Tissue_Label reference_tissue,
                           long seed, 
                           const char *rx_map_file, 
                           const char *tx_map_file) 
   : RF_Coil(seed, rx_map_file, tx_map_file){

   _noise_percentage     = percent;
   _reference_thickness  = reference_thickness;
   _reference_tissue     = reference_tissue;
}

//--------------------------------------------------------------------------
// Percent_Coil Destructor
//--------------------------------------------------------------------------

Percent_Coil::~Percent_Coil() {}

//--------------------------------------------------------------------------
// Percent_Coil::display_coil_info
// Writes information about the coil model to an output stream.
//--------------------------------------------------------------------------

void Percent_Coil::display_coil_info(ostream& stream) const {

   RF_Coil::display_coil_info(stream);

   if (_reference_thickness == 0) {
      stream << "Absolute percent noise coil." << endl;
      stream << "Percent noise:              " << _noise_percentage << "%" 
             << endl;
   } else {
      stream << "Reference percent noise coil." << endl;
      stream << "Percent noise:              " << _noise_percentage << "%"
             << endl;
      stream << "Reference slice thickness:  " << _reference_thickness 
             << " mm" << endl;
      stream << "Reference tissue label:     " << (int)_reference_tissue 
             << endl; 
   }
   stream << endl;

}

//--------------------------------------------------------------------------
// Percent_Coil::_compute_variance
// Computes the noise variance for the current pulse sequence.
//--------------------------------------------------------------------------

void Percent_Coil::_compute_variance(const Phantom &phantom){
   double ref_signal, std_dev;

   Pulse_Sequence &pseq = phantom.get_current_pulse_sequence();
   double slice_thickness = pseq.get_slice_thickness();

   // Compute the noise standard deviation as a percent of the
   // maximum signal acquired.
   // Compensate for varying slice thickness
   if (_reference_tissue == 0){
      ref_signal = phantom.get_max_mag_intensity();
   } else {
      ref_signal = phantom.get_mag_intensity(_reference_tissue);
   }

   if (_reference_thickness == 0){
      std_dev = (_noise_percentage/100.0)*ref_signal;
   } else {
      std_dev = (_noise_percentage/100.0)*ref_signal*
                (_reference_thickness/slice_thickness);
   }

   // Noise variance is the square of standard deviation.
   _noise_variance = SQR(std_dev);

}
