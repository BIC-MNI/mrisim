//===========================================================================
// DISCRETE_RF_PHANTOM.CXX
//
// R.Kwan
// April 25, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/discrete_rf_phantom.cxx,v 1.1 2003-05-30 16:43:10 bert Exp $
// $Log: discrete_rf_phantom.cxx,v $
// Revision 1.1  2003-05-30 16:43:10  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:09:42  rkwan
// Release 2.5
//
//===========================================================================

#include "discrete_rf_phantom.h"

//---------------------------------------------------------------------------
// Discrete_RF_Phantom constructor
//---------------------------------------------------------------------------

Discrete_RF_Phantom::Discrete_RF_Phantom(unsigned int n_tissue_classes,
                                         unsigned int n_flip_angles) :
   Phantom(n_tissue_classes),
   RF_Tissue_Phantom(n_tissue_classes, n_flip_angles),
   Discrete_Label_Phantom(n_tissue_classes) {

}

//---------------------------------------------------------------------------
// Discrete_RF_Phantom destructor
//---------------------------------------------------------------------------

Discrete_RF_Phantom::~Discrete_RF_Phantom() {

}

//---------------------------------------------------------------------------
// Discrete_RF_Phantom::use_rf_coil
// Associates an RF_Coil with the Phantom.
//---------------------------------------------------------------------------

int Discrete_RF_Phantom::use_rf_coil(const RF_Coil *rf_coil) {  

#ifdef DEBUG
   assert(rf_coil != NULL);
#endif

   int    status = TRUE; 
   double min, max;

   // Save a pointer to the RF Coil being used
   _rf_coil = (RF_Coil *)rf_coil;

   // Receive map checks and initialization

   if (rf_coil->uses_rx_map()) {
      if (!rf_coil->rx_map_is_same_size_as(_tissue_label_file)){
         cerr << "Receive map is not the same size as the Phantom." << endl;
         status = FALSE;
      }
   }

   // Transmit map checks and initialization

   if (rf_coil->uses_tx_map()) {
      if (!rf_coil->tx_map_is_same_size_as(_tissue_label_file)){
         cerr << "Transmit map is not the same size as the Phantom." << endl;
         status = FALSE;
      }
      rf_coil->get_tx_map_range(min, max);
      _setup_flip_errors(min, max);
   }

   return status;

}

//---------------------------------------------------------------------------
// Discrete_RF_Phantom::get_simulated_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Discrete_RF_Phantom::get_simulated_phantom_slice(int slice_num,
                                      Complex_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   // Clear the slice
   sim_slice.zeros();
 
   unsigned int m, n;
   Real_Slice rf_map(sim_slice.get_nrows(), sim_slice.get_ncols());
   MRI_Label  discrete_label(sim_slice.get_nrows(), sim_slice.get_ncols());

   // --- Load Labelled Phantom Data --- //

   _load_label_slice(slice_num, discrete_label);

   // --- Apply RF transmit inhomogeneity --- //

   if (uses_tx_map()) {

      // Load RF coil transmit map
      load_tx_map_slice(slice_num, rf_map);
   
      // Generate simulated slice with transmit inhomogeneity
      for (m=0; m<sim_slice.get_nrows(); m++){
         for (n=0; n<sim_slice.get_ncols(); n++){
   
            sim_slice.real(m,n) = RF_Tissue_Phantom::get_real_intensity(
                                    discrete_label(m,n), (double)rf_map(m,n) ); 
            sim_slice.imag(m,n) = RF_Tissue_Phantom::get_imag_intensity(
                                    discrete_label(m,n), (double)rf_map(m,n) );
         }
      }

   } else {
  
      // Generate simulated slice without transmit inhomogeneity
      for (m=0; m<sim_slice.get_nrows(); m++){
         for (n=0; n<sim_slice.get_ncols(); n++){

            sim_slice.real(m,n) = RF_Tissue_Phantom::get_real_intensity(
                                        discrete_label(m,n));
            sim_slice.imag(m,n) = RF_Tissue_Phantom::get_imag_intensity(
                                        discrete_label(m,n));

         }
      }  

   }

   
   // --- Apply RF receive inhomogeneity --- //

   // If a receive map has been specified, multiply the simulated slice
   // by the receive map, otherwise do nothing.

   if (uses_rx_map()) {

      load_rx_map_slice(slice_num, rf_map);
      sim_slice *= rf_map;
   
   }

}

//---------------------------------------------------------------------------
// Discrete_RF_Phantom::get_simulated_mag_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Discrete_RF_Phantom::get_simulated_mag_phantom_slice(int slice_num,
                                      Real_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   // Clear the slice
   sim_slice.zeros();
 
   unsigned int m, n;
   Real_Slice rf_map(sim_slice.get_nrows(), sim_slice.get_ncols());
   MRI_Label  discrete_label(sim_slice.get_nrows(), sim_slice.get_ncols());

   // --- Load Labelled Phantom Data --- //

   _load_label_slice(slice_num, discrete_label);

   // --- Apply RF transmit inhomogeneity --- //

   if (uses_tx_map()) {

      // Load RF coil transmit map
      load_tx_map_slice(slice_num, rf_map);
   
      // Generate simulated slice with transmit inhomogeneity
      for (m=0; m<sim_slice.get_nrows(); m++){
         for (n=0; n<sim_slice.get_ncols(); n++){
         
            sim_slice(m,n) = RF_Tissue_Phantom::get_mag_intensity(
                                   discrete_label(m,n), rf_map(m,n) );
         }
      }

   } else {
  
      // Generate simulated slice without transmit inhomogeneity
      for (m=0; m<sim_slice.get_nrows(); m++){
         for (n=0; n<sim_slice.get_ncols(); n++){

            sim_slice(m,n) = RF_Tissue_Phantom::get_mag_intensity(
                                   discrete_label(m,n));
         }
      }  

   }

   // --- Apply RF receive inhomogeneity --- //

   // If a receive map has been specified, multiply the simulated slice
   // by the receive map, otherwise do nothing.

   if (uses_rx_map()) {

      load_rx_map_slice(slice_num, rf_map);
      sim_slice *= rf_map;
   
   }

}

//---------------------------------------------------------------------------
// Discrete_RF_Phantom::get_simulated_real_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Discrete_RF_Phantom::get_simulated_real_phantom_slice(int slice_num,
                                      Real_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   // Clear the slice
   sim_slice.zeros();
 
   unsigned int m, n;
   Real_Slice rf_map(sim_slice.get_nrows(), sim_slice.get_ncols());
   MRI_Label  discrete_label(sim_slice.get_nrows(), sim_slice.get_ncols());

   // --- Load Labelled Phantom Data --- //

   _load_label_slice(slice_num, discrete_label);

   // --- Apply RF transmit inhomogeneity --- //

   if (uses_tx_map()) {

      // Load RF coil transmit map
      load_tx_map_slice(slice_num, rf_map);
   
      // Generate simulated slice with transmit inhomogeneity
      for (m=0; m<get_nrows(); m++){
         for (n=0; n<get_ncols(); n++){
   
            sim_slice(m,n) = RF_Tissue_Phantom::get_real_intensity(
                                   discrete_label(m,n), (double)rf_map(m,n) );
         }
      }

   } else {
  
      // Generate simulated slice without transmit inhomogeneity
      for (m=0; m<get_nrows(); m++){
         for (n=0; n<get_ncols(); n++){

            sim_slice(m,n) = RF_Tissue_Phantom::get_real_intensity(
                                   discrete_label(m,n));
         }
      }  

   }

   
   // --- Apply RF receive inhomogeneity --- //

   // If a receive map has been specified, multiply the simulated slice
   // by the receive map, otherwise do nothing.

   if (uses_rx_map()) {

      load_rx_map_slice(slice_num, rf_map);
      sim_slice *= rf_map;
   
   }

}
