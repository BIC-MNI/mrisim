//===========================================================================
// FUZZY_PHANTOM.CXX
//
// R.Kwan
// April 25, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/fuzzy_phantom.cxx,v 1.1 2003-05-30 16:43:10 bert Exp $
// $Log: fuzzy_phantom.cxx,v $
// Revision 1.1  2003-05-30 16:43:10  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  15:51:39  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:10:44  rkwan
// Release 2.5
//
//===========================================================================

#include "fuzzy_phantom.h"

//---------------------------------------------------------------------------
// Fuzzy_Phantom constructor
//---------------------------------------------------------------------------

Fuzzy_Phantom::Fuzzy_Phantom(unsigned int n_tissue_classes) :
   Phantom(n_tissue_classes), 
   Tissue_Phantom(n_tissue_classes),
   Fuzzy_Label_Phantom(n_tissue_classes) {

}

//---------------------------------------------------------------------------
// Fuzzy_Phantom destructor
//---------------------------------------------------------------------------

Fuzzy_Phantom::~Fuzzy_Phantom() {

}

//---------------------------------------------------------------------------
// Fuzzy_Phantom::get_simulated_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Fuzzy_Phantom::get_simulated_phantom_slice(int slice_num,
                                                Complex_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   // Create a fuzzy label slice of the same size
   Real_Slice fuzzy_label(sim_slice.get_nrows(), sim_slice.get_ncols());
   Real_Slice sum(sim_slice.get_nrows(), sim_slice.get_ncols());

   // Clear the slice
   sim_slice.zeros();
   sum.zeros();
 
   // Loop over each tissue type computing the weighting
   unsigned int itissue, m, n;
   Tissue_Label tissue_label;

   for (itissue=0; itissue<get_num_tissues(); itissue++){

      tissue_label = get_tissue_label(itissue);
      _load_label_slice(slice_num, tissue_label, fuzzy_label);

      //sim_slice.saxpy(get_real_intensity(tissue_label),
      //                get_imag_intensity(tissue_label), 
      //                fuzzy_label, sim_slice);

      for (m=0; m<sim_slice.get_nrows(); m++){
         for (n=0; n<sim_slice.get_ncols(); n++){
            sim_slice.real(m,n) += fuzzy_label(m,n) *
                     Tissue_Phantom::get_real_intensity(tissue_label);
            sim_slice.imag(m,n) += fuzzy_label(m,n) *
                     Tissue_Phantom::get_imag_intensity(tissue_label);
            sum(m,n) += fuzzy_label(m,n);
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

}

//---------------------------------------------------------------------------
// Fuzzy_Phantom::get_simulated_mag_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Fuzzy_Phantom::get_simulated_mag_phantom_slice(int slice_num,
                                                Real_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   // Create a fuzzy label slice of the same size
   Real_Slice fuzzy_label(sim_slice.get_nrows(), sim_slice.get_ncols());
   Real_Slice sum(sim_slice.get_nrows(), sim_slice.get_ncols());

   // Clear the slice
   sim_slice.zeros();
   sum.zeros();
 
   // Loop over each tissue type computing the weighting
   unsigned int itissue, m, n;
   Tissue_Label tissue_label;

   for (itissue=0; itissue<get_num_tissues(); itissue++){

      tissue_label = get_tissue_label(itissue);
      _load_label_slice(slice_num, tissue_label, fuzzy_label);

      for (m=0; m<sim_slice.get_nrows(); m++){
         for (n=0; n<sim_slice.get_ncols(); n++){
            sim_slice(m,n) += fuzzy_label(m,n) *
                              Tissue_Phantom::get_mag_intensity(tissue_label);
            sum(m,n) += fuzzy_label(m,n);
         }
      }

   }
 
   // --- Fuzzy volume normalization --- //
   sim_slice /= sum;

}

//---------------------------------------------------------------------------
// Fuzzy_Phantom::get_simulated_real_phantom_slice
// Generate a coloured phantom slice from labelled phantom data.
//---------------------------------------------------------------------------

void Fuzzy_Phantom::get_simulated_real_phantom_slice(int slice_num,
                                                Real_Slice& sim_slice) {

#ifdef DEBUG
   // Check that the slice is the right size
   assert(this->is_same_slice_size_as(sim_slice));
   assert(slice_num >= 0);
   assert(slice_num < get_nslices());
#endif

   // Create a fuzzy label slice of the same size
   Real_Slice fuzzy_label(sim_slice.get_nrows(), sim_slice.get_ncols());
   Real_Slice sum(sim_slice.get_nrows(), sim_slice.get_ncols());

   // Clear the slice
   sim_slice.zeros();
   sum.zeros();

   // Loop over each tissue type computing the weighting
   unsigned int itissue, m, n;
   Tissue_Label tissue_label;

   for (itissue=0; itissue<get_num_tissues(); itissue++){

      tissue_label = get_tissue_label(itissue);
      _load_label_slice(slice_num, tissue_label, fuzzy_label);

      //sim_slice.saxpy(get_real_intensity(tissue_label),
      //                fuzzy_label, sim_slice);

      for (m=0; m<sim_slice.get_nrows(); m++){
         for (n=0; n<sim_slice.get_ncols(); n++){
            sim_slice(m,n) += fuzzy_label(m,n) *
                              Tissue_Phantom::get_real_intensity(tissue_label);
            sum(m,n) += fuzzy_label(m,n);
         }
      }

   }

   // --- Fuzzy volume normalization --- //
   sim_slice /= sum;

}
