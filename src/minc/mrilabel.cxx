//===========================================================================
// MRILABEL.CXX
//
// R.Kwan
// April 19, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/mrilabel.cxx,v 1.1 2003-05-30 16:43:09 bert Exp $
// $Log: mrilabel.cxx,v $
// Revision 1.1  2003-05-30 16:43:09  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  15:48:16  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:23:09  rkwan
// Release 2.5
//
//
//===========================================================================

#include "mrilabel.h"
#include <math.h>

//---------------------------------------------------------------------------
// MRI_Label::get_mask(MRI_Byte_Matrix&)
// Returns a mask for a given label.
//---------------------------------------------------------------------------

void MRI_Label::get_mask(unsigned char label, MRI_Byte_Matrix& mask) const {
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(mask));
#endif

   for (n=0; n<len; n++){
      if (_matrix[n] == label){
         mask._matrix[n] = 1;
      } else {
         mask._matrix[n] = 0;
      }
   }
}

//---------------------------------------------------------------------------
// MRI_Label::get_mask(MRI_Short_Matrix&)
// Returns a mask for a given label.
//---------------------------------------------------------------------------

void MRI_Label::get_mask(unsigned char label, MRI_Short_Matrix& mask) const {
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(mask));
#endif

   for (n=0; n<len; n++){
      if (_matrix[n] == label){
         mask._matrix[n] = 1;
      } else {
         mask._matrix[n] = 0;
      }
   }
}

//---------------------------------------------------------------------------
// MRI_Label::get_mask(MRI_Float_Matrix&)
// Returns a mask for a given label.
//---------------------------------------------------------------------------

void MRI_Label::get_mask(unsigned char label, MRI_Float_Matrix& mask) const {
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(mask));
#endif

   for (n=0; n<len; n++){
      if (_matrix[n] == label){
         mask._matrix[n] = 1.0;
      } else {
         mask._matrix[n] = 0.0;
      }
   }
}

//---------------------------------------------------------------------------
// MRI_Label::get_mask(MRI_Double_Matrix&)
// Returns a mask for a given label.
//---------------------------------------------------------------------------

void MRI_Label::get_mask(unsigned char label, MRI_Double_Matrix& mask) const {
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(mask));
#endif

   for (n=0; n<len; n++){
      if (_matrix[n] == label){
         mask._matrix[n] = 1.0;
      } else {
         mask._matrix[n] = 0.0;
      }
   }
}

   
   
