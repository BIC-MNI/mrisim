#ifndef __IMAGE_SNR_COIL_H
#define __IMAGE_SNR_COIL_H

//==========================================================================
// IMAGE_SNR_COIL.H
// Image SNR coil class Image_SNR_Coil.
// Inherits from:  RF_Coil
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/image_snr_coil.h,v 1.1 2003-05-30 16:43:11 bert Exp $
 * $Log: image_snr_coil.h,v $
 * Revision 1.1  2003-05-30 16:43:11  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:11:47  rkwan
 * *** empty log message ***
 *
 * Revision 1.3  1995/12/22  20:13:57  rkwan
 * Update for percent coil and RF map features.
 *
 * Revision 1.2  1995/12/11  15:30:59  rkwan
 * Doc update
 *
 *========================================================================*/

#include "rf_coil.h"

//--------------------------------------------------------------------------
// Image_SNR_Coil Class
//--------------------------------------------------------------------------

class Image_SNR_Coil : public RF_Coil {
   public:
      Image_SNR_Coil(long seed = 0, 
                     const char *rx_map_file = NULL,
                     const char *tx_map_file = NULL);
      Image_SNR_Coil(double snr, long seed = 0, 
                     const char *rx_map_file = NULL,
                     const char *tx_map_file = NULL);

      virtual ~Image_SNR_Coil();

      // --- Random number generation --- //

      inline double get_image_snr(void) const;
      inline void   set_image_snr(double snr);

      // --- Convenience functions --- //

      virtual void display_coil_info(ostream& stream) const;
      virtual void _compute_variance(const Phantom &phantom);

   protected:

      // --- Internal data structures --- //
      double _image_snr;
};

//--------------------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Image_SNR_Coil::get_image_snr
// Returns the image SNR of the RF_Coil.
// Negative SNR indicates noiseless SNR.
//--------------------------------------------------------------------------

inline 
double Image_SNR_Coil::get_image_snr(void) const {
   return _image_snr;
}

//--------------------------------------------------------------------------
// Image_SNR_Coil::set_image_snr
// Sets the image SNR of the RF_Coil.
// Negative SNR indicates noiseless SNR.
//--------------------------------------------------------------------------

inline 
void   Image_SNR_Coil::set_image_snr(double snr) {
   _image_snr = snr;
}

#endif
