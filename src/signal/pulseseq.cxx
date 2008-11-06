//===========================================================================
// PULSESEQ.CXX
//
// R.Kwan
// August 1, 1995
//
// (C) Copyright 1995, 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/signal/pulseseq.cxx,v 1.2 2008-11-06 10:58:23 rotor Exp $
// $Log: pulseseq.cxx,v $
// Revision 1.2  2008-11-06 10:58:23  rotor
//  * fixed includes for iostream and friends
//  * updated for new release (1.0.2)
//
// Revision 1.1  2003/05/30 16:43:13  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  15:43:02  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:33:17  rkwan
// Release 2.5
//
//===========================================================================

#include "pulseseq.h"
#include <iostream>
#include <iomanip>
#include <string.h>

//---------------------------------------------------------------------------
// Pulse_Sequence constructors
//---------------------------------------------------------------------------

Pulse_Sequence::Pulse_Sequence(){

   _FOV[0]         = _FOV[1]         = _FOV[2]         = UNDEFINED;
   _matrix_size[0] = _matrix_size[1] = _matrix_size[2] = UNDEFINED;
   _offset[0]      = _offset[1]      = _offset[2]      = 0.0;

   _foldover_suppression   = NO;
   _foldover_direction     = COLUMN;
   _partial_fourier_method = NO_PARTIAL_FOURIER;
   _scan_percentage        = 1.0;
   _num_of_averages        = 1;
   _image_orientation      = TRANSVERSE;
   _image_type             = MOD_IMAGE;
   _scan_mode              = SCAN_MODE_MS;
   _water_fat_shift        = 2.0;

}

Pulse_Sequence::Pulse_Sequence(const Pulse_Sequence& p){
   int i;

   for(i=0; i<3; i++){
      _FOV[i]         = p._FOV[i];
      _matrix_size[i] = p._matrix_size[i];
      _offset[i]      = p._offset[i];
   }
   _foldover_suppression   = p._foldover_suppression;
   _foldover_direction     = p._foldover_direction;
   _partial_fourier_method = p._partial_fourier_method;
   _scan_percentage        = p._scan_percentage;
   _num_of_averages        = p._num_of_averages;
   _image_orientation      = p._image_orientation;
   _image_type             = p._image_type;
   _scan_mode              = p._scan_mode;
   _water_fat_shift        = p._water_fat_shift;

}

//---------------------------------------------------------------------------
// Pulse_Sequence destructor
//---------------------------------------------------------------------------

Pulse_Sequence::~Pulse_Sequence(){
   // Do nothing
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_voxel_step
// Returns the step size of the given dimension.
//--------------------------------------------------------------------------

double Pulse_Sequence::get_voxel_step(int i) const {

#ifdef DEBUG
   assert(i >= 0 && i < 3);
   assert(_matrix_size[i] != 0);
#endif

   double step;
   if (_FOV[i] == UNDEFINED || _matrix_size[i] == UNDEFINED) {
      step = _step[i];
   } else if (i == SLICE) {
      step = (_FOV[SLICE] - _slice_thickness) / (_matrix_size[SLICE] - 1);
   } else {
      step = _FOV[i] / (double)_matrix_size[i];
   }
   return step;
}

void Pulse_Sequence::get_voxel_step(double step[]) const {
   step[0] = get_voxel_step(0);
   step[1] = get_voxel_step(1);
   step[2] = get_voxel_step(2);
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_voxel_step
// Sets the voxel step if FOV and matrix size are not used.
//--------------------------------------------------------------------------

void Pulse_Sequence::set_voxel_step(double step[]) { 
   if (_FOV[0] == UNDEFINED || _matrix_size[0] == UNDEFINED) {
      _step[0] = step[0];
   }
   if (_FOV[1] == UNDEFINED || _matrix_size[1] == UNDEFINED) {
      _step[1] = step[1];
   }
   if (_FOV[2] == UNDEFINED || _matrix_size[2] == UNDEFINED) {
      _step[2] = step[2];
   }
}
   
//--------------------------------------------------------------------------
// Pulse_Sequence::get_axis
// Returns coordinate -> axis lookup depending on image orientation.
//--------------------------------------------------------------------------

Axis Pulse_Sequence::get_axis(int coordinate) const {
   static Axis coord_transformation[3][3] = {
      {Z_AXIS, Y_AXIS, X_AXIS},  // TRANSVERSE
      {X_AXIS, Z_AXIS, Y_AXIS},  // SAGITTAL
      {Y_AXIS, Z_AXIS, X_AXIS}}; // CORONAL

   return coord_transformation[_image_type][coordinate];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_partial_fourier_name
// Returns the name of the given partial fourier method.
//--------------------------------------------------------------------------

const char *Pulse_Sequence::get_partial_fourier_name(
    Partial_Fourier_Type method) const {

   static const char *method_name[] = {"NO_PARTIAL_FOURIER", 
                                      "HALF_FOURIER",
                                      "PARTIAL_ECHO",
                                      "PARTIAL_MATRIX",
                                      "RECTANGULAR_FOV" };

   return method_name[method];

}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_orientation_name
// Returns the name of the given image orientation.
//--------------------------------------------------------------------------

const char *Pulse_Sequence::get_orientation_name(
   Image_Orientation orientation) const {

   static const char *orientation_name[] = {"SAME", "TRANSVERSE",
                                           "SAGITTAL", "CORONAL" };

   return orientation_name[orientation];

}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_image_type_name
// Returns the name of the given image type.
//--------------------------------------------------------------------------

const char *Pulse_Sequence::get_image_type_name(Image_Type image_type) const {

   static const char *image_type_name[] = {"NONE","M","R","P","I"};
   return image_type_name[image_type];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_scan_mode_name
// Returns the name of the given scan mode.
//--------------------------------------------------------------------------

const char *Pulse_Sequence::get_scan_mode_name(Scan_Mode scan_mode) const {
   static const char *scan_mode_name[] = {"MS", "3D", "2D" };
   return scan_mode_name[scan_mode];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::display_info
// Output information about the pulse sequence.
//--------------------------------------------------------------------------

void Pulse_Sequence::display_info(ostream& stream) const {

   const char *partial_fourier_method = 
      get_partial_fourier_name(_partial_fourier_method);
   const char *image_orientation = 
      get_orientation_name(_image_orientation);
   const char *image_type =
      get_image_type_name(_image_type);
   const char *scan_mode = 
      get_scan_mode_name(_scan_mode);

   stream << "Pulse Sequence Information:" << endl;
   stream << "---------------------------" << endl << endl;

   stream << "SLICE:  ";
   stream << "FOV: " << setw(6) << get_fov(SLICE) << " mm ";
   stream << "Matrix Size: " << setw(4) << get_matrix_size(SLICE) << " ";
   stream << "Voxel Step: " << setw(3) << get_voxel_step(SLICE) << " ";
   stream << "Offcentre: " << setw(3) << get_voxel_offset(SLICE) << " mm" 
          << endl;

   stream << "ROW:    ";
   stream << "FOV: " << setw(6) << get_fov(ROW) << " mm ";
   stream << "Matrix Size: " << setw(4) << get_matrix_size(ROW) << " ";
   stream << "Voxel Step: " << setw(3) << get_voxel_step(ROW) << " ";
   stream << "Offcentre: " << setw(3) << get_voxel_offset(ROW) << " mm" 
          << endl;

   stream << "COLUMN: ";
   stream << "FOV: " << setw(6) << get_fov(COLUMN) << " mm ";
   stream << "Matrix Size: " << setw(4) << get_matrix_size(COLUMN) << " ";
   stream << "Voxel Step: " << setw(3) << get_voxel_step(COLUMN) << " ";
   stream << "Offcentre: " << setw(3) << get_voxel_offset(COLUMN) << " mm" 
          << endl << endl;

   stream << "Foldover_suppression = " 
          << (_foldover_suppression ? "YES" : "NO") << endl;
   stream << "Foldover direction   = "
          << (_foldover_direction == COLUMN ? "X" : "Y") << endl;
   stream << "Partial Fourier      = " << partial_fourier_method << endl;
   stream << "Scan Percentage      = " << get_scan_percentage() << endl;
   stream << "Number of Averages   = " << get_num_of_averages() << endl;
   stream << "Image Orientation    = " << image_orientation << endl;
   stream << "Image Display Type   = " << image_type << endl;
   stream << "Scan Mode            = " << scan_mode << endl;
   stream << "Water Fat Shift      = " << get_water_fat_shift() 
          << " pixels " << endl;
   stream << "Sampling Period      = " << get_sampling_period()
          << " us" << endl;
   stream << endl;

}
