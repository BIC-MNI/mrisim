#ifndef __MRIIMAGE_H
#define __MRIIMAGE_H

//===========================================================================
// MRIIMAGE.H
// Image class.
// Inherits from:  MRI_Short_Matrix
// Base class to:  
// 
// R.Kwan
// September 16, 1995
//
// (C) Copyright 1995 by R.Kwan
//===========================================================================

/*===========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/mriimage.h,v 1.1 2003-05-30 16:43:09 bert Exp $
 * $Log: mriimage.h,v $
 * Revision 1.1  2003-05-30 16:43:09  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:47:35  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:22:43  rkwan
 * Release 2.5
 *
 * Revision 1.4  1995/12/22  20:17:58  rkwan
 * Update for percent coil and RF map features.
 *
 * Revision 1.3  1995/12/11  15:15:05  rkwan
 * Fix RCS header bug.
 *
 *=========================================================================*/

#include "mrimatrix.h"

//---------------------------------------------------------------------------
// MRI_Image class
// Specialization of MRI_Short_Matrix to store slices of a MINC volume.
//---------------------------------------------------------------------------

class MRI_Image : public MRI_Short_Matrix {
   public:
      MRI_Image(unsigned int nrows, unsigned int ncols);

      MRI_Image(const MRI_Byte_Matrix& mat);
      MRI_Image(const MRI_Float_Matrix& mat);
      MRI_Image(const MRI_Double_Matrix& mat);

      virtual ~MRI_Image();

      // --- VIO_Real range access functions --- //

      inline double get_real_minimum(void) const;
      inline double get_real_maximum(void) const;
             void   set_real_minimum(double min);
             void   set_real_maximum(double max);

      // --- Image pixel manipulation --- //

      inline short  convert_value_to_voxel(double value) const;
      inline double convert_voxel_to_value(short voxel) const;

      inline short  get_voxel(unsigned int row, unsigned int col) const;
      inline double get_value(unsigned int row, unsigned int col) const;
      inline void   set_voxel(unsigned int row, unsigned int col, 
                              short voxel);
      inline void   set_value(unsigned int row, unsigned int col,
                              double value);

      // --- Image Scaling --- //
 
      MRI_Image& operator+=(short a);
      MRI_Image& operator*=(double a);
      MRI_Image& operator*=(const MRI_Float_Matrix &mat);

      MRI_Image& scale(double scale);
      MRI_Image& offset(double offset);

   private:

      MRI_Image();

      // --- Internal data members --- //
      double _real_minimum;
      double _real_maximum;
      short  _voxel_minimum;
      short  _voxel_maximum;

      double _voxel_to_real_scale;
      double _real_to_voxel_scale;

      // --- Internal member functions --- //
      inline void _update_scaling_factors(void);

};

//---------------------------------------------------------------------------
// MRI_Image::get_real_minimum
// Returns the minimum voxel value in real values.
//---------------------------------------------------------------------------

inline
double MRI_Image::get_real_minimum(void) const {
   return _real_minimum;
}

//---------------------------------------------------------------------------
// MRI_Image::get_real_maximum
// Returns the maximum voxel value in real values.
//---------------------------------------------------------------------------

inline
double MRI_Image::get_real_maximum(void) const {
   return _real_maximum;
}

//---------------------------------------------------------------------------
// MRI_Image::convert_value_to_voxel
// Converts a real value to a voxel value.
//---------------------------------------------------------------------------

inline
short MRI_Image::convert_value_to_voxel(double value) const {

#ifdef DEBUG
/*
   //assert(value >= _real_minimum - 1E-5);
   if (value < _real_minimum) {
      cout << value << " " << _real_minimum << endl;
      abort();
   }   
   //assert(value <= _real_maximum + 1E-5);
   if (value > _real_maximum) {
      cout << value << " " << _real_maximum << endl;
      abort();
   }   
*/
#endif

   short voxel;
   if (value > _real_maximum) {
      voxel = _voxel_maximum;
   } else if (value < _real_minimum) {
      voxel = _voxel_minimum;
   } else {
      voxel = (short) rint( _real_to_voxel_scale * (value - _real_minimum) +
                        _voxel_minimum );
   }

#ifdef DEBUG
   //assert(voxel >= _voxel_minimum);
   if (voxel < _voxel_minimum){
      cout << "min: " << value << " " << _real_minimum << " " << _real_maximum
           << " " << voxel << " " << _voxel_minimum << endl;
      abort();
   }
   //assert(voxel <= _voxel_maximum);
   if (voxel > _voxel_maximum){
      cout << "max: " << value << " " << _real_minimum << " " << _real_maximum
           << " " << voxel << " " << _voxel_maximum << endl;
      abort();
   }
#endif

   return voxel;
}

//---------------------------------------------------------------------------
// MRI_Image::convert_voxel_to_value
// Converts a voxel value to a real value.
//---------------------------------------------------------------------------

inline 
double MRI_Image::convert_voxel_to_value(short voxel) const {

#ifdef DEBUG
   assert(voxel >= _voxel_minimum);
   assert(voxel <= _voxel_maximum);
#endif

   double value;
   value = _voxel_to_real_scale*(voxel-_voxel_minimum) + _real_minimum; 

#ifdef DEBUG
   assert(value >= _real_minimum - 1E-5);
   assert(value <= _real_maximum + 1E-5);
#endif

   return value;
}

//---------------------------------------------------------------------------
// MRI_Image::get_voxel
// Returns the actual short image voxel value.
//---------------------------------------------------------------------------

inline 
short MRI_Image::get_voxel(unsigned int row, unsigned int col) const {
   return _matrix[row*_ncols+col];
}

//---------------------------------------------------------------------------
// MRI_Image::get_value
// Returns the scaled real double image voxel value.
//---------------------------------------------------------------------------

inline 
double MRI_Image::get_value(unsigned int row, unsigned int col) const {
   return convert_voxel_to_value(_matrix[row*_ncols+col]);
}

//---------------------------------------------------------------------------
// MRI_Image::set_voxel
// Stores an actual short image voxel in the image. 
//---------------------------------------------------------------------------

inline 
void MRI_Image::set_voxel(unsigned int row, unsigned int col, short voxel) {
   _matrix[row*_ncols+col] = voxel;
}

//---------------------------------------------------------------------------
// MRI_Image::set_value
// Scales and stores a real value in the image.
//---------------------------------------------------------------------------

inline 
void MRI_Image::set_value(unsigned int row, unsigned int col, double value) {
   _matrix[row*_ncols+col] = convert_value_to_voxel(value);
}

//---------------------------------------------------------------------------
// MRI_Image::_update_scaling_factors
// Recomputes the voxel->real and real->voxel scaling factors.
//---------------------------------------------------------------------------

inline 
void MRI_Image::_update_scaling_factors(void) {
   _voxel_to_real_scale = (_real_maximum - _real_minimum) /
                          (double)(_voxel_maximum - _voxel_minimum);
   _real_to_voxel_scale = 1.0 / _voxel_to_real_scale;
}

#endif

