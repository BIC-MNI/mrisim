#ifndef __DISCRETE_RF_PHANTOM_H
#define __DISCRETE_RF_PHANTOM_H

//==========================================================================
// DISCRETE_RF_PHANTOM.H
// Discrete_RF_Phantom class.
// Inherits from:  Discrete_Label_Phantom, RF_Tissue_Phantom
// Base class to:  
//
// R. Kwan
// (C) Copyright 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/discrete_rf_phantom.h,v 1.1 2003-05-30 16:43:10 bert Exp $
 * $Log: discrete_rf_phantom.h,v $
 * Revision 1.1  2003-05-30 16:43:10  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:09:33  rkwan
 * Release 2.5
 *
 *========================================================================*/

#include "discrete_label_phantom.h"
#include "rf_tissue_phantom.h"

//---------------------------------------------------------------------------
// Discrete_RF_Phantom class
//---------------------------------------------------------------------------

class Discrete_RF_Phantom : public Discrete_Label_Phantom, 
                            public RF_Tissue_Phantom {
   public:
      Discrete_RF_Phantom(unsigned int n_tissue_classes,
                          unsigned int n_flip_angles);

      virtual ~Discrete_RF_Phantom();

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

   private:
      
      // --- Internal member functions --- //
      inline void load_safe_rx_map_slice(int slice_num, 
                                         Real_Slice& rx_slice);
      inline void load_safe_tx_map_slice(int slice_num, 
                                         Real_Slice& tx_slice);

};

//---------------------------------------------------------------------------
// Inline member functions
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Discrete_RF_Phantom::has_same_dimension_info_as_rx 
// Returns TRUE if the receive map has the same dimension info as the
// Phantom.
//---------------------------------------------------------------------------

inline 
int Discrete_RF_Phantom::has_same_dimension_info_as_rx(
      const RF_Coil &rf_coil) const{
   return rf_coil.rx_map_has_same_dimension_info_as(_tissue_label_file);
}

//---------------------------------------------------------------------------
// Discrete_RF_Phantom::has_same_dimension_info_as_tx 
// Returns TRUE if the transmit map has the same dimension info as the
// Phantom.
//---------------------------------------------------------------------------

inline 
int Discrete_RF_Phantom::has_same_dimension_info_as_tx(
      const RF_Coil &rf_coil) const{
   return rf_coil.tx_map_has_same_dimension_info_as(_tissue_label_file);
}

//---------------------------------------------------------------------------
// Discrete_RF_Phantom::has_same_dimension_names_as_rx 
// Returns TRUE if the receive map has the same dimension names as the
// Phantom.
//---------------------------------------------------------------------------

inline 
int Discrete_RF_Phantom::has_same_dimension_names_as_rx(
      const RF_Coil &rf_coil) const {
   return rf_coil.rx_map_has_same_dimension_names_as(_tissue_label_file);
}

//---------------------------------------------------------------------------
// Discrete_RF_Phantom::has_same_dimension_names_as_tx 
// Returns TRUE if the transmit map has the same dimension names as the
// Phantom.
//---------------------------------------------------------------------------

inline 
int Discrete_RF_Phantom::has_same_dimension_names_as_tx(
      const RF_Coil &rf_coil) const {
   return rf_coil.tx_map_has_same_dimension_names_as(_tissue_label_file);
}

//---------------------------------------------------------------------------
// Discrete_RF_Phantom::load_safe_rx_map_slice
// Loads a receive coil inhomogeneity slice.
// Checks if the slice number is valid, if it is outside the valid range
// sets the slice to zero.
//---------------------------------------------------------------------------

inline
void Discrete_RF_Phantom::load_safe_rx_map_slice(int slice_num,
                                                 Real_Slice& rx_slice){

   if ((slice_num >= 0) && (slice_num < get_nslices())) {
      load_rx_map_slice(slice_num, rx_slice);
   } else {
      rx_slice.zeros();
   }

}

//---------------------------------------------------------------------------
// Discrete_RF_Phantom::load_safe_tx_map_slice
// Loads a transmit coil inhomogeneity slice.
// Checks if the slice number is valid, if it is outside the valid range
// sets the slice to zero.
//---------------------------------------------------------------------------

inline
void Discrete_RF_Phantom::load_safe_tx_map_slice(int slice_num,
                                                 Real_Slice& tx_slice){

   if ((slice_num >= 0) && (slice_num < get_nslices())) {
      load_tx_map_slice(slice_num, tx_slice);
   } else {
      tx_slice.zeros();
   }

}

#endif

