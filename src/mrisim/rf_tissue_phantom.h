#ifndef __RF_TISSUE_PHANTOM_H
#define __RF_TISSUE_PHANTOM_H

//===========================================================================
// RF_TISSUE_PHANTOM.H
// RF_Tissue_Phantom class.
// Inherits from:  Phantom
// Base class to:  Discrete_RF_Phantom, Fuzzy_RF_Phantom
//
// R.Kwan
// April 29, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

/*===========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/rf_tissue_phantom.h,v 1.1 2003-05-30 16:43:12 bert Exp $
 * $Log: rf_tissue_phantom.h,v $
 * Revision 1.1  2003-05-30 16:43:12  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:19:18  rkwan
 * Release 2.5
 *
 *=========================================================================*/

#include "phantom.h"
#include "rf_coil.h"
#include <signal/quickseq.h>
#include <signal/customseq.h>
#include <signal/tissue.h>

//---------------------------------------------------------------------------
// RF_Tissue_Phantom class
// Implementation base class that stores information about tissue types
// in the Phantom, and implements interface to Pulse_Sequence for signal
// production simulation with RF inhomogeneities.
//---------------------------------------------------------------------------

class RF_Tissue_Phantom : virtual public Phantom {
   public:
      RF_Tissue_Phantom(unsigned int n_tissue_classes,
                        unsigned int n_flip_angles);

      virtual ~RF_Tissue_Phantom();

      // --- Tissue information --- //
      void display_tissue_info(ostream& stream) const;

      // --- Pulse Sequence simulation --- //
      void apply_pulse_sequence(Quick_Sequence *pseq);
      void find_steady_state(Custom_Sequence *pseq);

      // --- Access functions --- //
      inline int uses_rx_map(void) const; 
      inline int uses_tx_map(void) const;

      inline void load_rx_map_slice(int slice_num, Real_Slice& rf_slice);
      inline void load_tx_map_slice(int slice_num, Real_Slice& tx_slice);

      inline double get_real_intensity(Tissue_Label label, 
                                       double flip_error) const;
      inline double get_imag_intensity(Tissue_Label label, 
                                       double flip_error) const;
      inline double get_mag_intensity(Tissue_Label label, 
                                      double flip_error) const;

      inline double get_real_intensity(Tissue_Label label) const;
      inline double get_imag_intensity(Tissue_Label label) const;
      inline double get_mag_intensity(Tissue_Label label) const;

   protected:

      void get_signal_map(unsigned int tissue_index,
                          const Real_Slice& tx_map,
                          Complex_Slice& signal_map);

      void _setup_flip_errors(double min, double max);

      // --- RF Coil interface --- //
      RF_Coil       *_rf_coil;

   private:

      // --- Internal member functions --- //
      double _interp_real_intensity(unsigned int index,
                                double flip_error) const;
      double _interp_imag_intensity(unsigned int index,
                                double flip_error) const;
      double _interp_mag_intensity(unsigned int index,
                               double flip_error) const;

      inline double& _lookup_real_intensity(unsigned int tissue_index,
                                 unsigned int flip_index) const;
      inline double& _lookup_imag_intensity(unsigned int tissue_index,
                                 unsigned int flip_index) const;

      // --- Internal data structures --- //
      unsigned int  _n_flip_angles;
      double        *_flip_error;     // Flip angle errors (stored by index)
      double        _flip_error_step;

      double        *_real_intensity; // Computed signal intensities
      double        *_imag_intensity; // (stored by index)

      double        *_no_error_real;
      double        *_no_error_imag;

};

//---------------------------------------------------------------------------
// Inline member functions
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::_lookup_real_intensity
// Lookup computed real tissue intensity.
//---------------------------------------------------------------------------

inline
double& RF_Tissue_Phantom::_lookup_real_intensity(unsigned int tissue_index,
                                 unsigned int flip_index) const {

   return _real_intensity[tissue_index*_n_flip_angles+flip_index];
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::_lookup_imag_intensity
// Lookup computed imag tissue intensity.
//---------------------------------------------------------------------------

inline
double& RF_Tissue_Phantom::_lookup_imag_intensity(unsigned int tissue_index,
                                 unsigned int flip_index) const {

   return _imag_intensity[tissue_index*_n_flip_angles+flip_index];
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::get_real_intensity
// Lookup computed real tissue intensity.
//---------------------------------------------------------------------------

inline
double RF_Tissue_Phantom::get_real_intensity(Tissue_Label label,
                                             double flip_error) const {
   return _interp_real_intensity(get_tissue_index(label), flip_error);
}

inline
double RF_Tissue_Phantom::get_real_intensity(Tissue_Label label) const {
   return _no_error_real[get_tissue_index(label)];
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::get_imag_intensity
// Lookup computed imag tissue intensity.
//---------------------------------------------------------------------------

inline
double RF_Tissue_Phantom::get_imag_intensity(Tissue_Label label,
                                             double flip_error) const {
   return _interp_imag_intensity(get_tissue_index(label), flip_error);
}

inline
double RF_Tissue_Phantom::get_imag_intensity(Tissue_Label label) const {
   return _no_error_imag[get_tissue_index(label)];
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::get_mag_intensity
// Lookup computed magnitude tissue intensity.
//---------------------------------------------------------------------------

inline
double RF_Tissue_Phantom::get_mag_intensity(Tissue_Label label,
                                            double flip_error) const {
   return _interp_mag_intensity(get_tissue_index(label), flip_error);
}

inline
double RF_Tissue_Phantom::get_mag_intensity(Tissue_Label label) const {
   return hypot( _no_error_real[get_tissue_index(label)],
                 _no_error_imag[get_tissue_index(label)] );
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::uses_rx_map
// Returns TRUE if a receive coil inhomogeneity map was specified.
//---------------------------------------------------------------------------

inline
int RF_Tissue_Phantom::uses_rx_map(void) const {
   int status;
   if (_rf_coil != NULL) {
      status = _rf_coil->uses_rx_map();
   } else {
      status = FALSE;
   } 
   return status;
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::uses_tx_map
// Returns TRUE if a transmit coil inhomogeneity map was specified.
//---------------------------------------------------------------------------

inline
int RF_Tissue_Phantom::uses_tx_map(void) const {
   int status;
   if (_rf_coil != NULL) {
      status = _rf_coil->uses_tx_map();
   } else {
      status = FALSE;
   }
   return status;
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::load_rx_map_slice
// Loads a slice of the receive coil map.
//---------------------------------------------------------------------------

inline 
void RF_Tissue_Phantom::load_rx_map_slice(int slice_num, Real_Slice& rx_slice) {
#ifdef DEBUG
   assert(uses_rx_map());
#endif
   _rf_coil->load_rx_map_slice(slice_num, rx_slice);
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::load_tx_map_slice
// Loads a slice of the transmit coil map.
//---------------------------------------------------------------------------

inline 
void RF_Tissue_Phantom::load_tx_map_slice(int slice_num, Real_Slice& tx_slice) {
#ifdef DEBUG
   assert(uses_tx_map());
#endif
   _rf_coil->load_tx_map_slice(slice_num, tx_slice);
}

#endif
