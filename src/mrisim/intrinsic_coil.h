#ifndef __INTRINSIC_COIL_H
#define __INTRINSIC_COIL_H

//==========================================================================
// INTRINSIC_COIL.H
// Intrinsic SNR Coil class Intrinsic_SNR_Coil.
// Inherits from:  RF_Coil
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/intrinsic_coil.h,v 1.1 2003-05-30 16:43:11 bert Exp $
 * $Log: intrinsic_coil.h,v $
 * Revision 1.1  2003-05-30 16:43:11  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:12:40  rkwan
 * Release 2.5
 *
 * Revision 1.3  1995/12/22  20:15:01  rkwan
 * Update for percent coil and RF map features.
 *
 * Revision 1.2  1995/12/11  15:31:08  rkwan
 * Doc update.
 *
 *========================================================================*/

#include "rf_coil.h"

//--------------------------------------------------------------------------
// Intrinsic_SNR_Coil Class
//--------------------------------------------------------------------------

class Intrinsic_SNR_Coil : public RF_Coil {
   public:
      Intrinsic_SNR_Coil(long seed = 0, 
                         const char *rx_map_file = NULL,
                         const char *tx_map_file = NULL);
      Intrinsic_SNR_Coil(double ISNR, long seed = 0, 
                         const char *rx_map_file = NULL,
                         const char *tx_map_file = NULL);

      virtual ~Intrinsic_SNR_Coil();

      // --- Random number generation --- //
      inline double get_intrinsic_snr(void) const;
      inline void   set_intrinsic_snr(double ISNR);

      // --- Convenience functions --- //

      virtual void display_coil_info(ostream& stream) const;
      virtual void _compute_variance(const Phantom &phantom);

   protected:

      // --- Internal data structures --- //
      double _intrinsic_snr;
};

//--------------------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Intrinsic_SNR_Coil::get_intrinsic_snr
// Returns the intrinsic SNR of the RF_Coil
// Negative intrinsic SNR indicates no noise.
//--------------------------------------------------------------------------

inline
double Intrinsic_SNR_Coil::get_intrinsic_snr(void) const {
   return _intrinsic_snr;
}

//--------------------------------------------------------------------------
// Intrinsic_SNR_Coil::set_intrinsic_snr
// Sets the intrinsic SNR of the RF_Coil.
// Negative intrinsic SNR indicates no noise.
//--------------------------------------------------------------------------

inline
void Intrinsic_SNR_Coil::set_intrinsic_snr(double ISNR){
   // Set the intrinsic SNR (Hz^(1/2)/mm^3).
   _intrinsic_snr = ISNR;
}

#endif
