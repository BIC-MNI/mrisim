//==========================================================================
// INTRINSIC_COIL.CXX
// Intrinsic SNR Coil class Intrinsic_SNR_Coil.
// Inherits from:  RF_Coil
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/intrinsic_coil.cxx,v 1.1 2003-05-30 16:43:11 bert Exp $
// $Log: intrinsic_coil.cxx,v $
// Revision 1.1  2003-05-30 16:43:11  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:12:45  rkwan
// Release 2.5
//
// Revision 1.3  1995/12/22  20:15:07  rkwan
// Update for percent coil and RF map features.
//
// Revision 1.2  1995/12/11  15:19:06  rkwan
// *** empty log message ***
//
//==========================================================================

#include "intrinsic_coil.h"

//--------------------------------------------------------------------------
// Intrinsic_SNR_Coil Class
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Intrinsic_SNR_Coil Constructors
//--------------------------------------------------------------------------

Intrinsic_SNR_Coil::Intrinsic_SNR_Coil(long seed, const char *rx_map_file,
                                       const char *tx_map_file) 
   : RF_Coil(seed, rx_map_file, tx_map_file) {
   _intrinsic_snr = -1;
}

Intrinsic_SNR_Coil::Intrinsic_SNR_Coil(double ISNR, long seed, 
                                       const char *rx_map_file, 
                                       const char *tx_map_file) 
   : RF_Coil(seed, rx_map_file, tx_map_file) {
   _intrinsic_snr = ISNR;
}

//--------------------------------------------------------------------------
// Intrinsic_SNR_Coil Destructor
//--------------------------------------------------------------------------

Intrinsic_SNR_Coil::~Intrinsic_SNR_Coil() {}

//--------------------------------------------------------------------------
// Intrinsic_SNR_Coil::display_coil_info
// Writes information about a coil model to an output stream.
//--------------------------------------------------------------------------

void Intrinsic_SNR_Coil::display_coil_info(ostream& stream) const {

   RF_Coil::display_coil_info(stream);

   stream << "Intrinsic SNR RF Coil." << endl;
   stream << "Intrinsic SNR:              " << _intrinsic_snr << endl;
   stream << endl;
}

//--------------------------------------------------------------------------
// Intrinsic_SNR_Coil::_compute_variance
// Computes the noise variance for the current pulse sequence.
// Must be called before any noise generation.
//--------------------------------------------------------------------------

void Intrinsic_SNR_Coil::_compute_variance(const Phantom &phantom){
   double voxel_volume, partial_fourier_factor;
   Time_us sampling_period;

   if (_intrinsic_snr < 0){

      // A negative intrinsic SNR indicates a noiseless condition.
      _noise_variance = 0;

   } else {

      // Compute noise variance based on pulse sequence parameters.

      Pulse_Sequence &pseq = phantom.get_current_pulse_sequence();
      voxel_volume = pseq.get_voxel_volume();
      sampling_period = pseq.get_sampling_period();
     
      switch(pseq.get_partial_fourier_method()){
         case PARTIAL_MATRIX:
            partial_fourier_factor = 1/(pseq.get_scan_percentage());
            break;
         case HALF_FOURIER:
         case PARTIAL_ECHO:
         case RECTANGULAR_FOV:
            partial_fourier_factor = pseq.get_scan_percentage();
            break;
         case NO_PARTIAL_FOURIER:
         default:
            partial_fourier_factor = 1.0;
            break;
      }

      if (pseq.uses_foldover_suppression()){
         partial_fourier_factor *= 2.0; 
      }
 
      if (pseq.get_scan_mode() == SCAN_MODE_3D){
         partial_fourier_factor *= pseq.get_matrix_size(SLICE);
      }

      _noise_variance = 1.0/
                        ((partial_fourier_factor)*
                         SQR(voxel_volume)*
                         SQR(_intrinsic_snr)*
                         sampling_period*1E-6*
                         pseq.get_num_of_averages());

   }
}
