#ifndef __PERCENT_COIL_H
#define __PERCENT_COIL_H

//==========================================================================
// PERCENT_COIL.H
// Percent SNR Coil class.
// Inherits from:  RF_Coil
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/percent_coil.h,v 1.1 2003-05-30 16:43:11 bert Exp $
 * $Log: percent_coil.h,v $
 * Revision 1.1  2003-05-30 16:43:11  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:17:32  rkwan
 * Release 2.5
 *
 * Revision 1.3  1995/12/22  20:15:14  rkwan
 * Update for percent coil and RF map features.
 *
 * Revision 1.2  1995/12/11  15:31:13  rkwan
 * Doc update.
 *
 *========================================================================*/

#include "rf_coil.h"

//--------------------------------------------------------------------------
// Percent_Coil Class
//--------------------------------------------------------------------------

class Percent_Coil : public RF_Coil {
   public:
      Percent_Coil(long         seed = 0, 
                   const char   *tx_map_file = NULL,
                   const char   *rx_map_file = NULL);

      Percent_Coil(double       percent, 
                   double       reference_thickness = 0,
                   Tissue_Label reference_tissue    = 0, 
                   long         seed                = 0, 
                   const char   *rx_map_file        = NULL, 
                   const char   *tx_map_file        = NULL);

      virtual ~Percent_Coil();

      // --- Random number generation --- //

      inline double       get_noise_percentage(void) const;
      inline double       get_reference_thickness(void) const;
      inline Tissue_Label get_reference_tissue_label(void) const;

      void set_noise_percentage(double percent);
      void set_reference_thickness(double reference_thickness);
      void set_reference_tissue_label(Tissue_Label label);

      // --- Convenience functions --- //

      virtual void display_coil_info(ostream& stream) const;
      virtual void _compute_variance(const Phantom &phantom);

   protected:

      // --- Internal data structures --- //
      double       _noise_percentage;
      double       _reference_thickness;
      Tissue_Label _reference_tissue;
};

//--------------------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Percent_Coil::get_noise_percentage
// Returns the noise percentage (0..100) of the RF_Coil
//--------------------------------------------------------------------------

inline
double Percent_Coil::get_noise_percentage(void) const {
   return _noise_percentage;
}

//--------------------------------------------------------------------------
// Percent_Coil::get_reference_thickness
// Returns the reference noise slice thickness (in mm)
//--------------------------------------------------------------------------

inline
double Percent_Coil::get_reference_thickness(void) const {
   return _reference_thickness;
}

//--------------------------------------------------------------------------
// Percent_Coil::get_reference_tissue
// Returns the reference noise tissue label.
//--------------------------------------------------------------------------

inline
Tissue_Label Percent_Coil::get_reference_tissue_label(void) const {
   return _reference_tissue;
}

//--------------------------------------------------------------------------
// Percent_Coil::set_noise_percentage
// Sets the noise percentage (0..100) of the RF_Coil
//--------------------------------------------------------------------------

inline
void Percent_Coil::set_noise_percentage(double percent){

#ifdef DEBUG
   // Ensure noise percentage is in the valid range.
   assert(percent >= 0.0);
   assert(percent <= 100.0);
#endif
   _noise_percentage = percent;
}

//--------------------------------------------------------------------------
// Percent_Coil::set_reference_thickness
// Sets the reference noise slice thickness (in mm).
//--------------------------------------------------------------------------

inline
void Percent_Coil::set_reference_thickness(double reference_thickness){
   _reference_thickness = reference_thickness;
}

//--------------------------------------------------------------------------
// Percent_Coil::set_reference_tissue_label
// Sets the reference noise tissue label.
//--------------------------------------------------------------------------

inline
void Percent_Coil::set_reference_tissue_label(Tissue_Label label) {
   _reference_tissue = label;
}

#endif
