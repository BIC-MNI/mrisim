#ifndef __TISSUE_PHANTOM_H
#define __TISSUE_PHANTOM_H

//===========================================================================
// TISSUE_PHANTOM.H
// Tissue_Phantom class.
// Inherits from:  Phantom
// Base class to:  Discrete_Phantom, Fuzzy_Phantom
//
// R.Kwan
// April 29, 1996
//
// (C) Copyright 1996 by R. Kwan
//===========================================================================

/*===========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/tissue_phantom.h,v 1.1 2003-05-30 16:43:12 bert Exp $
 * $Log: tissue_phantom.h,v $
 * Revision 1.1  2003-05-30 16:43:12  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:19:36  rkwan
 * Release 2.5
 *
 *=========================================================================*/

#include "phantom.h"
#include <signal/quickseq.h>
#include <signal/customseq.h>
#include <signal/tissue.h>

//---------------------------------------------------------------------------
// Tissue_Phantom class
// Implementation base class that stores information about tissue types
// in the Phantom, and implements interface to Pulse_Sequence for signal
// production simulation.
//---------------------------------------------------------------------------

class Tissue_Phantom : virtual public Phantom {
   public:
      Tissue_Phantom(unsigned int n_tissue_classes);

      virtual ~Tissue_Phantom();

      // --- Tissue information --- //
      void display_tissue_info(ostream& stream) const;

      // --- Pulse Sequence simulation --- //
      void apply_pulse_sequence(Quick_Sequence *pseq);
      void find_steady_state(Custom_Sequence *pseq);

      // --- Access functions --- //
      inline double get_real_intensity(Tissue_Label label) const;
      inline double get_imag_intensity(Tissue_Label label) const;
      inline double get_mag_intensity(Tissue_Label label) const;

   private:

      inline double _lookup_real_intensity(unsigned int index) const;
      inline double _lookup_imag_intensity(unsigned int index) const;
      inline double _lookup_mag_intensity(unsigned int index) const;

      // --- Internal data structures --- //
      double        *_real_intensity; // Computed signal intensities
      double        *_imag_intensity; // (stored by index)

};

//---------------------------------------------------------------------------
// Inline member functions
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Tissue_Phantom::get_real_intensity
// Lookup the real intensity for the given tissue label.
//---------------------------------------------------------------------------

inline
double Tissue_Phantom::get_real_intensity(Tissue_Label label) const {
   return _lookup_real_intensity(get_tissue_index(label));
}

//---------------------------------------------------------------------------
// Tissue_Phantom::get_imag_intensity
// Lookup the imag intensity for the given tissue label.
//---------------------------------------------------------------------------

inline
double Tissue_Phantom::get_imag_intensity(Tissue_Label label) const {
   return _lookup_imag_intensity(get_tissue_index(label));
}

//---------------------------------------------------------------------------
// Tissue_Phantom::get_mag_intensity
// Lookup the magnitude intensity for the given tissue label.
//---------------------------------------------------------------------------

inline
double Tissue_Phantom::get_mag_intensity(Tissue_Label label) const {
   return _lookup_mag_intensity(get_tissue_index(label));
}

//---------------------------------------------------------------------------
// Private member functions
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Tissue_Phantom::_lookup_real_intensity
// Returns the simulated intensity for a tissue index.
//---------------------------------------------------------------------------

inline
double Tissue_Phantom::_lookup_real_intensity(unsigned int index) const {
   return _real_intensity[index];
}

//---------------------------------------------------------------------------
// Tissue_Phantom::_lookup_imag_intensity
// Returns the simulated intensity for a tissue index.
//---------------------------------------------------------------------------

inline
double Tissue_Phantom::_lookup_imag_intensity(unsigned int index) const {
   return _imag_intensity[index];
}

//---------------------------------------------------------------------------
// Tissue_Phantom::_lookup_mag_intensity
// Returns the simulated intensity for a tissue index.
//---------------------------------------------------------------------------

inline
double Tissue_Phantom::_lookup_mag_intensity(unsigned int index) const {
   return hypot(_real_intensity[index], _imag_intensity[index]);
}

#endif


