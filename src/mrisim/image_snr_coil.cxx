//==========================================================================
// IMAGE_SNR_COIL.CXX
// Image SNR coil class Image_SNR_Coil.
// Inherits from:  RF_Coil
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/image_snr_coil.cxx,v 1.1 2003-05-30 16:43:11 bert Exp $
// $Log: image_snr_coil.cxx,v $
// Revision 1.1  2003-05-30 16:43:11  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:12:14  rkwan
// Release 2.5
//
// Revision 1.3  1995/12/22  20:14:50  rkwan
// Update for percent coil and RF map features.
//
// Revision 1.2  1995/12/11  15:18:42  rkwan
// *** empty log message ***
//
//==========================================================================

#include "image_snr_coil.h"

#ifdef DEBUG
#include <assert.h>
#endif

//--------------------------------------------------------------------------
// Image_SNR_Coil Class
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Image_SNR_Coil Constructors
//--------------------------------------------------------------------------

Image_SNR_Coil::Image_SNR_Coil(long seed, const char *rx_map_file,
                               const char *tx_map_file) 
   : RF_Coil(seed, rx_map_file, tx_map_file) {
   _image_snr = -1;
}

Image_SNR_Coil::Image_SNR_Coil(double snr, long seed, 
                               const char *rx_map_file, 
                               const char *tx_map_file) 
   : RF_Coil(seed, rx_map_file, tx_map_file) {
   _image_snr = snr;
}

//--------------------------------------------------------------------------
// Image_SNR_Coil Destructor
//--------------------------------------------------------------------------

Image_SNR_Coil::~Image_SNR_Coil() {}

//--------------------------------------------------------------------------
// Image_SNR_Coil::display_coil_info
// Writes information about a coil model to an output stream.
//--------------------------------------------------------------------------

void Image_SNR_Coil::display_coil_info(ostream& stream) const {

   RF_Coil::display_coil_info(stream);

   stream << "Image SNR Coil." << endl;
   stream << "Image SNR:                  " << _image_snr << endl;
   stream << endl;

}

//--------------------------------------------------------------------------
// Image_SNR_Coil::_compute_variance
// Computes the noise variance for the current pulse sequence.
//--------------------------------------------------------------------------

void Image_SNR_Coil::_compute_variance(const Phantom &phantom){
   double max_signal, std_dev;

   // Compute the noise standard deviation as a percent of the
   // maximum signal acquired.

   max_signal = phantom.get_max_mag_intensity();
   std_dev = max_signal/(_image_snr*0.66);

   // Noise variance is the square of standard deviation.
   _noise_variance = SQR(std_dev);
}
