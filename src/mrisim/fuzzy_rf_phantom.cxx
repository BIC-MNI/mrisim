//===========================================================================
// FUZZY_RF_PHANTOM.CXX
//
// R.Kwan
// April 25, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/fuzzy_rf_phantom.cxx,v 1.1 2003-05-30 16:43:10 bert Exp $
// $Log: fuzzy_rf_phantom.cxx,v $
// Revision 1.1  2003-05-30 16:43:10  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  15:56:00  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:11:06  rkwan
// Release 2.5
//
//===========================================================================

#include "fuzzy_rf_phantom.h"

//---------------------------------------------------------------------------
// Fuzzy_RF_Phantom constructor
//---------------------------------------------------------------------------

Fuzzy_RF_Phantom::Fuzzy_RF_Phantom(unsigned int n_tissue_classes,
                                         unsigned int n_flip_angles) :
   Phantom(n_tissue_classes),
   RF_Tissue_Phantom(n_tissue_classes, n_flip_angles),
   Fuzzy_Label_Phantom(n_tissue_classes) {

}

//---------------------------------------------------------------------------
// Fuzzy_RF_Phantom destructor
//---------------------------------------------------------------------------

Fuzzy_RF_Phantom::~Fuzzy_RF_Phantom() {

}

//---------------------------------------------------------------------------
// Fuzzy_RF_Phantom::use_rf_coil
// Associates an RF_Coil with the Phantom.
//---------------------------------------------------------------------------

int Fuzzy_RF_Phantom::use_rf_coil(const RF_Coil *rf_coil) {       

#ifdef DEBUG
   assert(rf_coil != NULL);
#endif

   int    status = TRUE; 
   double min, max;

   // Save a pointer to the RF Coil being used
   _rf_coil = (RF_Coil *)rf_coil;

   // Receive map checks and initialization

   if (rf_coil->uses_rx_map()) {
      if (!rf_coil->rx_map_is_same_size_as(_tissue_label_file[0])){
         cerr << "Receive map is not the same size as the Phantom." << endl;
         status = FALSE;
      }
   }

   // Transmit map checks and initialization

   if (rf_coil->uses_tx_map()) {
      if (!rf_coil->tx_map_is_same_size_as(_tissue_label_file[0])){
         cerr << "Transmit map is not the same size as the Phantom." << endl;
         status = FALSE;
      }
      rf_coil->get_tx_map_range(min, max);
      _setup_flip_errors(min, max);
   }

   return status;

}

//---------------------------------------------------------------------------
// Fuzzy_RF_Phantom::get_simulated_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Fuzzy_RF_Phantom::get_simulated_phantom_slice(int slice_num,
                                      Complex_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   Real_Slice rf_map(sim_slice.get_nrows(), sim_slice.get_ncols());
   Real_Slice fuzzy_label(sim_slice.get_nrows(), sim_slice.get_ncols());
   Real_Slice sum(sim_slice.get_nrows(), sim_slice.get_ncols());

   // Clear slice
   sim_slice.zeros();
   sum.zeros();

   // --- Apply RF transmit inhomogeneity --- //

   unsigned int itissue, m, n;
   Tissue_Label tissue_label;

   if (uses_tx_map()) {

      // Load RF coil transmit map
      load_tx_map_slice(slice_num, rf_map);
 
      // Generate simulated slice with transmit inhomogeneity 
      // Loop over each tissue type computing signal contribution
      for (itissue=0; itissue<get_num_tissues(); itissue++){

         tissue_label = get_tissue_label(itissue);
         _load_label_slice(slice_num, tissue_label, fuzzy_label);

         for (m=0; m<sim_slice.get_nrows(); m++){
            for (n=0; n<sim_slice.get_ncols(); n++){

               sim_slice.real(m,n) += fuzzy_label(m,n) *
                        RF_Tissue_Phantom::get_real_intensity(tissue_label,
                                                         (double)rf_map(m,n));
               sim_slice.imag(m,n) += fuzzy_label(m,n) *
                        RF_Tissue_Phantom::get_imag_intensity(tissue_label,
                                                         (double)rf_map(m,n));
               sum(m,n) += fuzzy_label(m,n);
            }
         }

      }
   
   } else {

      // Generate simulated slice without transmit inhomogeneity
      // Loop over each tissue type computing signal contribution
      for (itissue=0; itissue<get_num_tissues(); itissue++){
 
         tissue_label = get_tissue_label(itissue);
         _load_label_slice(slice_num, tissue_label, fuzzy_label);

         for (m=0; m<sim_slice.get_nrows(); m++){
            for (n=0; n<sim_slice.get_ncols(); n++){
   
               sim_slice.real(m,n) += fuzzy_label(m,n) *
                        RF_Tissue_Phantom::get_real_intensity(tissue_label);
               sim_slice.imag(m,n) += fuzzy_label(m,n) *
                        RF_Tissue_Phantom::get_imag_intensity(tissue_label);
               sum(m,n) += fuzzy_label(m,n);

            }
         }

      }

   }

   // --- Fuzzy volume normalization --- //

   for (m=0; m<sim_slice.get_nrows(); m++){
      for (n=0; n<sim_slice.get_ncols(); n++){
         sim_slice.real(m,n) /= sum(m,n);
         sim_slice.imag(m,n) /= sum(m,n);
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
// Fuzzy_RF_Phantom::get_simulated_mag_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Fuzzy_RF_Phantom::get_simulated_mag_phantom_slice(int slice_num,
                                      Real_Slice& sim_slice) {


#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   Real_Slice rf_map(sim_slice.get_nrows(), sim_slice.get_ncols());
   Real_Slice fuzzy_label(sim_slice.get_nrows(), sim_slice.get_ncols());
   Real_Slice sum(sim_slice.get_nrows(), sim_slice.get_ncols());

   // Clear the slice
   sim_slice.zeros();
   sum.zeros();

   // --- Apply RF transmit inhomogeneity --- //

   unsigned int itissue, m, n;
   Tissue_Label tissue_label;

   if (uses_tx_map()) {

      // Load RF coil transmit map
      load_tx_map_slice(slice_num, rf_map);

      // Generate simulated slice with transmit inhomogeneity
      // Loop over each tissue type computing signal contribution
      for (itissue=0; itissue<get_num_tissues(); itissue++){

         tissue_label = get_tissue_label(itissue);
         _load_label_slice(slice_num, tissue_label, fuzzy_label);

         for (m=0; m<sim_slice.get_nrows(); m++){
            for (n=0; n<sim_slice.get_ncols(); n++){

               sim_slice(m,n) += fuzzy_label(m,n) *
                        RF_Tissue_Phantom::get_mag_intensity(tissue_label,
                                                         (double)rf_map(m,n));
               sum(m,n) += fuzzy_label(m,n);
            }
         }

      }

   } else {

      // Generate simulated slice without transmit inhomogeneity
      // Loop over each tissue type computing signal contribution
      for (itissue=0; itissue<get_num_tissues(); itissue++){

         tissue_label = get_tissue_label(itissue);
         _load_label_slice(slice_num, tissue_label, fuzzy_label);

         for (m=0; m<sim_slice.get_nrows(); m++){
            for (n=0; n<sim_slice.get_ncols(); n++){

               sim_slice(m,n) += fuzzy_label(m,n) *
                        RF_Tissue_Phantom::get_mag_intensity(tissue_label);
               sum(m,n) += fuzzy_label(m,n);
            }
         }

      }

   }

   // --- Fuzzy volume normalization --- // 
   sim_slice /= sum;

   // --- Apply RF receive inhomogeneity --- //

   // If a receive map has been specified, multiply the simulated slice
   // by the receive map, otherwise do nothing.

   if (uses_rx_map()) {

      load_rx_map_slice(slice_num, rf_map);
      sim_slice *= rf_map;

   }

}

//---------------------------------------------------------------------------
// Fuzzy_RF_Phantom::get_simulated_real_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Fuzzy_RF_Phantom::get_simulated_real_phantom_slice(int slice_num,
                                      Real_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   Real_Slice rf_map(sim_slice.get_nrows(), sim_slice.get_ncols());
   Real_Slice fuzzy_label(sim_slice.get_nrows(), sim_slice.get_ncols());
   Real_Slice sum(sim_slice.get_nrows(), sim_slice.get_ncols());

   // Clear the slice
   sim_slice.zeros();
   sum.zeros();

   // --- Apply RF transmit inhomogeneity --- //

   unsigned int itissue, m, n;
   Tissue_Label tissue_label;

   if (uses_tx_map()) {

      // Load RF coil transmit map
      load_tx_map_slice(slice_num, rf_map);

      // Generate simulated slice with transmit inhomogeneity
      // Loop over each tissue type computing signal contribution
      for (itissue=0; itissue<get_num_tissues(); itissue++){

         tissue_label = get_tissue_label(itissue);
         _load_label_slice(slice_num, tissue_label, fuzzy_label);

         for (m=0; m<sim_slice.get_nrows(); m++){
            for (n=0; n<sim_slice.get_ncols(); n++){

               sim_slice(m,n) += fuzzy_label(m,n) *
                        RF_Tissue_Phantom::get_real_intensity(tissue_label,
                                                         (double)rf_map(m,n));
               sum(m,n) += fuzzy_label(m,n);
            }
         }

      }

   } else {

      // Generate simulated slice without transmit inhomogeneity
      // Loop over each tissue type computing signal contribution
      for (itissue=0; itissue<get_num_tissues(); itissue++){

         tissue_label = get_tissue_label(itissue);
         _load_label_slice(slice_num, tissue_label, fuzzy_label);

         for (m=0; m<sim_slice.get_nrows(); m++){
            for (n=0; n<sim_slice.get_ncols(); n++){

               sim_slice(m,n) += fuzzy_label(m,n) *
                        RF_Tissue_Phantom::get_real_intensity(tissue_label);
               sum(m,n) += fuzzy_label(m,n);
            }
         }

      }

   }

   // --- Fuzzy volume normalization --- // 
   sim_slice /= sum;

   // --- Apply RF receive inhomogeneity --- //

   // If a receive map has been specified, multiply the simulated slice
   // by the receive map, otherwise do nothing.

   if (uses_rx_map()) {

      load_rx_map_slice(slice_num, rf_map);
      sim_slice *= rf_map;

   }

}
