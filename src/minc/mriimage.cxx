//===========================================================================
// MRIIMAGE.CXX
//
// R.Kwan
// October 16, 1995
//
// (C) Copyright 1995 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/mriimage.cxx,v 1.1 2003-05-30 16:43:09 bert Exp $
// $Log: mriimage.cxx,v $
// Revision 1.1  2003-05-30 16:43:09  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.6  1996/05/29  18:53:53  rkwan
// *** empty log message ***
//
// Revision 2.5.1.1  1996/05/29  18:44:36  rkwan
// Fix set_real_minimum/set_real_maximum
//
// Revision 2.5  1996/05/29  16:22:51  rkwan
// Release 2.5
//
// Revision 1.3  1995/12/22  20:17:52  rkwan
// Update for percent coil and RF map features.
//
// Revision 1.2  1995/12/11  14:16:47  rkwan
// Updated for fast_iso_model.
//
//===========================================================================

#include "mriimage.h"
#include <math.h>

//---------------------------------------------------------------------------
// MRI_Image constructors
//---------------------------------------------------------------------------

MRI_Image::MRI_Image() : MRI_Short_Matrix() {
   _voxel_minimum = 0;
   _voxel_maximum = 4095;
   _real_minimum  = 0.0;
   _real_maximum  = 1.0;

   _update_scaling_factors();
}

MRI_Image::MRI_Image(unsigned int nrows, unsigned int ncols) 
   : MRI_Short_Matrix(nrows, ncols) {
   _voxel_minimum = 0;
   _voxel_maximum = 4095;
   _real_minimum  = 0.0;
   _real_maximum  = 1.0;

   _update_scaling_factors();
}

//---------------------------------------------------------------------------
// MRI_Image destructor 
//---------------------------------------------------------------------------

MRI_Image::~MRI_Image() {}

//---------------------------------------------------------------------------
// MRI_Image type cast constructor
//---------------------------------------------------------------------------

MRI_Image::MRI_Image(const MRI_Byte_Matrix& mat) :
   MRI_Short_Matrix(mat) {

   _voxel_minimum = 0;
   _voxel_maximum = 4095;
   _real_minimum  = (double)mat.minimum();
   _real_maximum  = (double)mat.maximum();

   _update_scaling_factors();
}

MRI_Image::MRI_Image(const MRI_Float_Matrix& mat) :
   MRI_Short_Matrix(mat.get_nrows(), mat.get_ncols()) {
   _voxel_minimum = 0;
   _voxel_maximum = 4095;
   _real_minimum  = mat.minimum();
   _real_maximum  = mat.maximum();

   _update_scaling_factors();

   unsigned int n;
   unsigned int len = get_nelements();

   for(n=0; n<len; n++){
      _matrix[n] = convert_value_to_voxel(mat._matrix[n]);
   }

}

MRI_Image::MRI_Image(const MRI_Double_Matrix& mat) :
   MRI_Short_Matrix(mat.get_nrows(), mat.get_ncols()) {

   _voxel_minimum = 0;
   _voxel_maximum = 4095;
   _real_minimum  = mat.minimum();
   _real_maximum  = mat.maximum();

   _update_scaling_factors();

   unsigned int n;
   unsigned int len = get_nelements();

   for(n=0; n<len; n++){
      _matrix[n] = convert_value_to_voxel(mat._matrix[n]);
   }

}

//---------------------------------------------------------------------------
// MRI_Image::set_real_minimum
// Sets the real maximum value for the image.
//---------------------------------------------------------------------------

void MRI_Image::set_real_minimum(double min) {
   _real_minimum = min;
   _update_scaling_factors();
}

//---------------------------------------------------------------------------
// MRI_Image:set_real_maximum
// Sets the real maximum value for the image.
//---------------------------------------------------------------------------

void MRI_Image::set_real_maximum(double max) {
   _real_maximum = max;
   _update_scaling_factors();
}

//---------------------------------------------------------------------------
// MRI_Image::operator+=
// Adds an offset to all elements in the image.
//---------------------------------------------------------------------------

MRI_Image& MRI_Image::operator+=(short a){
   unsigned int n;
   unsigned int len = this->get_nelements();

   for (n=0; n<len; n++){
      _matrix[n] += a;
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Image::operator*=(const double)
// Scales all elements in the image.
//---------------------------------------------------------------------------

MRI_Image& MRI_Image::operator*=(double a){
   unsigned int n;
   unsigned int len = this->get_nelements();

   for (n=0; n<len; n++){
      _matrix[n] = (short) rint( _matrix[n]*a );
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Image::operator*=(const mriFloatMatrix&)
// Piecewise multiplication of all elements in the matrix.
//---------------------------------------------------------------------------

MRI_Image& MRI_Image::operator*=(const MRI_Float_Matrix& mat){

#ifdef DEBUG
   assert(this->is_same_size_as(mat));
#endif

   unsigned int m, n;

   for (m=0; m<get_nrows(); m++){
      for (n=0; n<get_ncols(); n++){
         (*this)(m,n) = (short) rint ( (*this)(m,n) * mat(m,n) );
      }
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Image::scale
// Scales all elements in the image.
//---------------------------------------------------------------------------

MRI_Image& MRI_Image::scale(double scale){
   _real_minimum *= scale;
   _real_maximum *= scale;

   _update_scaling_factors();

   return *this;
}

//---------------------------------------------------------------------------
// MRI_Image::offset
// Offsets all elements in the image.
//---------------------------------------------------------------------------

MRI_Image& MRI_Image::offset(double offset) {
   _real_minimum += offset;
   _real_maximum += offset;

   _update_scaling_factors();

   return *this;
}

