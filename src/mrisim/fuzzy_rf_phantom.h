#ifndef __FUZZY_RF_PHANTOM_H
#define __FUZZY_RF_PHANTOM_H

//==========================================================================
// FUZZY_RF_PHANTOM.H
// Fuzzy_RF_Phantom class.
// Inherits from:  Fuzzy_Label_Phantom, RF_Tissue_Phantom
// Base class to:  
//
// R. Kwan
// (C) Copyright 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/fuzzy_rf_phantom.h,v 1.1 2003-05-30 16:43:10 bert Exp $
 * $Log: fuzzy_rf_phantom.h,v $
 * Revision 1.1  2003-05-30 16:43:10  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:52:11  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:10:58  rkwan
 * Release 2.5
 *
 *========================================================================*/

#include "fuzzy_label_phantom.h"
#include "rf_tissue_phantom.h"

//---------------------------------------------------------------------------
// Fuzzy_RF_Phantom class
//---------------------------------------------------------------------------

class Fuzzy_RF_Phantom : public Fuzzy_Label_Phantom, 
                         public RF_Tissue_Phantom {
   public:
      Fuzzy_RF_Phantom(unsigned int n_tissue_classes,
                       unsigned int n_flip_angles);

      virtual ~Fuzzy_RF_Phantom();

      // --- RF Coil Information --- //
      int use_rf_coil(const RF_Coil *rf_coil);

      inline int has_same_dimension_info_as_rx(const RF_Coil &rf_coil) const;
      inline int has_same_dimension_info_as_tx(const RF_Coil &rf_coil) const;
      inline int has_same_dimension_names_as_rx(const RF_Coil &rf_coil) const;
      inline int has_same_dimension_names_as_tx(const RF_Coil &rf_coil) const;

      // --- Image Formation --- // 
      void get_simulated_phantom_slice(int slice_num,
                                       Complex_Slice& sim_slice);
      void get_simulated_mag_phantom_slice(int slice_num,
                                       Real_Slice& sim_slice);
      void get_simulated_real_phantom_slice(int slice_num,
                                       Real_Slice& sim_slice);

};

//---------------------------------------------------------------------------
// Inline member functions
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Fuzzy_RF_Phantom::has_same_dimension_info_as_rx
// Returns TRUE if the receive map has the same dimension info as the
// Phantom.
//---------------------------------------------------------------------------

inline
int Fuzzy_RF_Phantom::has_same_dimension_info_as_rx(
      const RF_Coil &rf_coil) const {
   return rf_coil.rx_map_has_same_dimension_info_as(_tissue_label_file[0]);
}

//---------------------------------------------------------------------------
// Fuzzy_RF_Phantom::has_same_dimension_info_as_tx
// Returns TRUE if the transmit map has the same dimension info as the
// Phantom.
//---------------------------------------------------------------------------

inline
int Fuzzy_RF_Phantom::has_same_dimension_info_as_tx(
      const RF_Coil &rf_coil) const {
   return rf_coil.tx_map_has_same_dimension_info_as(_tissue_label_file[0]);
}

//---------------------------------------------------------------------------
// Fuzzy_RF_Phantom::has_same_dimension_names_as_rx
// Returns TRUE if the receive map has the same dimension names as the
// Phantom.
//---------------------------------------------------------------------------

inline
int Fuzzy_RF_Phantom::has_same_dimension_names_as_rx(
      const RF_Coil &rf_coil) const {
   return rf_coil.rx_map_has_same_dimension_names_as(_tissue_label_file[0]);
}

//---------------------------------------------------------------------------
// Fuzzy_RF_Phantom::has_same_dimension_names_as_tx
// Returns TRUE if the transmit map has the same dimension names as the
// Phantom.
//---------------------------------------------------------------------------

inline
int Fuzzy_RF_Phantom::has_same_dimension_names_as_tx(
      const RF_Coil &rf_coil) const {
   return rf_coil.tx_map_has_same_dimension_names_as(_tissue_label_file[0]);
}

#endif

