//===========================================================================
// DISCRETE_PHANTOM.CXX
//
// R.Kwan
// April 25, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/discrete_phantom.cxx,v 1.1 2003-05-30 16:43:10 bert Exp $
// $Log: discrete_phantom.cxx,v $
// Revision 1.1  2003-05-30 16:43:10  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:09:18  rkwan
// Release 2.5
//
//===========================================================================

#include "discrete_phantom.h"

//---------------------------------------------------------------------------
// Discrete_Phantom constructor
//---------------------------------------------------------------------------

Discrete_Phantom::Discrete_Phantom(unsigned int n_tissue_classes) :
   Phantom(n_tissue_classes),
   Tissue_Phantom(n_tissue_classes), 
   Discrete_Label_Phantom(n_tissue_classes) {

}

//---------------------------------------------------------------------------
// Discrete_Phantom destructor
//---------------------------------------------------------------------------

Discrete_Phantom::~Discrete_Phantom() {

}

//---------------------------------------------------------------------------
// Discrete_Phantom::get_simulated_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Discrete_Phantom::get_simulated_phantom_slice(int slice_num,
                                         Complex_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   // Create a discrete label slice for the phantom
   MRI_Label discrete_label(sim_slice.get_nrows(), sim_slice.get_ncols());
   _load_label_slice(slice_num, discrete_label);

   unsigned int m, n;
   for (m=0; m<sim_slice.get_nrows(); m++){
      for (n=0; n<sim_slice.get_ncols(); n++){
         sim_slice.real(m,n) = Tissue_Phantom::get_real_intensity(
                                     discrete_label(m,n));
         sim_slice.imag(m,n) = Tissue_Phantom::get_imag_intensity(
                                     discrete_label(m,n));
      }
   }

/*
   // --- Second Algorithm:  matrix operations --- //

   // Clear the slice
   sim_slice.zeros();

   // Create a discrete label slice and mask of the same size
   MRI_Label       discrete_label(get_nrows(), get_ncols());
   MRI_Byte_Matrix tissue_mask(get_nrows(), get_ncols());
   _load_label_slice(slice_num, discrete_label);

   // Loop over each tissue type computing the weighting
   unsigned int itissue;
   for (itissue=0; itissue<get_num_tissues(); itissue++){
      // Generate a mask for the tissue type
      discrete_label.get_mask(get_tissue_label(itissue), tissue_mask);
      sim_slice.saxpy(Tissue_Phantom::get_real_intensity(
                            get_tissue_label(itissue)),
                      Tissue_Phantom::get_imag_intensity(
                            get_tissue_label(itissue)),
                      tissue_mask, sim_slice);
   }

*/

}

//---------------------------------------------------------------------------
// Discrete_Phantom::get_simulated_mag_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Discrete_Phantom::get_simulated_mag_phantom_slice(int slice_num,
                                         Real_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   // Create a discrete label slice for the phantom
   MRI_Label discrete_label(sim_slice.get_nrows(), sim_slice.get_ncols());
   _load_label_slice(slice_num, discrete_label);

   unsigned int m, n;
   for (m=0; m<sim_slice.get_nrows(); m++){
      for (n=0; n<sim_slice.get_ncols(); n++){
         sim_slice(m,n) = Tissue_Phantom::get_mag_intensity(
                                discrete_label(m,n));
      }
   }
}

//---------------------------------------------------------------------------
// Discrete_Phantom::get_simulated_real_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Discrete_Phantom::get_simulated_real_phantom_slice(int slice_num,
                                         Real_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   // Create a discrete label slice for the phantom
   MRI_Label discrete_label(sim_slice.get_nrows(), sim_slice.get_ncols());
   _load_label_slice(slice_num, discrete_label);

   unsigned int m, n;
   for (m=0; m<sim_slice.get_nrows(); m++){
      for (n=0; n<sim_slice.get_ncols(); n++){
         sim_slice(m,n) = Tissue_Phantom::get_real_intensity(
                                discrete_label(m,n));
      }
   }
}
