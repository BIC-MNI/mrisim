#ifndef __PULSESEQ_H
#define __PULSESEQ_H

//==========================================================================
// PULSESEQ.H
//
// Base class for Pulse Sequence classes
// Inherits from:
// Base class to:  Quick_Sequence, Custom_Sequence
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/signal/pulseseq.h,v 1.1 2003-05-30 16:43:13 bert Exp $
 * $Log: pulseseq.h,v $
 * Revision 1.1  2003-05-30 16:43:13  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:42:27  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:33:12  rkwan
 * Release 2.5
 *
 *========================================================================*/

#include <mrisim/mrisim.h>
#include <stdio.h>
#include "vector.h"

//--------------------------------------------------------------------------
// Specialized Pulse Sequence parameter types
//--------------------------------------------------------------------------

enum Image_Type {NO_IMAGE = 0, MOD_IMAGE = 1, REAL_IMAGE = 2, 
                 PHASE_IMAGE = 3, IMAG_IMAGE = 4};

enum Image_Orientation {SAME = 0, TRANSVERSE = 1, SAGITTAL = 2, CORONAL = 3};

enum Scan_Mode {SCAN_MODE_MS = 0, SCAN_MODE_3D = 1, SCAN_MODE_2D = 2};

enum Partial_Fourier_Type {NO_PARTIAL_FOURIER = 0, HALF_FOURIER = 1,
                           PARTIAL_ECHO = 2, PARTIAL_MATRIX = 3,
                           RECTANGULAR_FOV = 4};

//--------------------------------------------------------------------------
// Image Option specification
//--------------------------------------------------------------------------

#ifndef YES
#define YES TRUE
#define NO  FALSE
#endif

#ifndef TRUE
#define TRUE  TRUE
#define FALSE FALSE
#endif

#ifndef UNDEFINED
#define UNDEFINED -1
#endif

//--------------------------------------------------------------------------
// Image axis specification
//--------------------------------------------------------------------------

#ifndef SLICE
#define SLICE  0
#define ROW    1
#define COLUMN 2
#endif

#define SLICE_SELECT SLICE
#define PHASE_ENCODE ROW
#define READOUT      COLUMN

//--------------------------------------------------------------------------
// Pulse_Sequence class
// Base class for Pulse Sequence classes.
//--------------------------------------------------------------------------

class Pulse_Sequence {
   public:
      Pulse_Sequence();                          // constructors
      Pulse_Sequence(const Pulse_Sequence& p);

      virtual ~Pulse_Sequence();                 // destructor

      // --- Voxel Geometry Options --- //

      inline int     get_matrix_size(int i) const;
      inline void    get_matrix_size(int matrix_size[]) const;
      inline int     get_reconstruction_size(int i) const;
      inline int     get_nslices(void) const;
      inline int     get_nrows(void) const;
      inline int     get_ncols(void) const;
      inline double  get_fov(int i) const;
      inline void    get_fov(double fov[]) const;
      inline double  get_voxel_offset(int i) const;
      inline double  get_voxel_volume(void) const;
             double  get_voxel_step(int i) const;
             void    get_voxel_step(double step[]) const;
      inline double  get_slice_thickness(void) const;

      inline void    set_fov(double fov[]);
      inline void    set_matrix_size(int matrix_size[]);
      inline void    set_voxel_offset(double offset[]);
             void    set_voxel_step(double step[]);
      inline void    set_slice_thickness(double slice_thickness);
      inline void    set_reconstruction_size(int reconstruction_size[]);

             Axis    get_axis(int coordinate) const;

      // --- Partial Fourier Options --- //
      
      inline int uses_foldover_suppression(void) const;
      inline int uses_partial_fourier(void) const;

      inline Partial_Fourier_Type get_partial_fourier_method(void) const;
      inline double               get_scan_percentage(void) const;
      inline int                  get_foldover_direction(void) const;

      inline void set_foldover_suppression(int on_off);
      inline void set_foldover_direction(int direction);
      inline void set_partial_fourier_method(Partial_Fourier_Type method, 
                                             double scan_percentage);
     
      // --- Image Acquisition Options --- //
 
      inline int               get_num_of_averages(void) const;
      inline Image_Orientation get_image_orientation(void) const;
      inline Image_Type        get_image_type(void) const;
      inline Scan_Mode         get_scan_mode(void) const;

      inline void set_num_of_averages(int num_of_averages);
      inline void set_image_orientation(Image_Orientation orientation);
      inline void set_image_type(Image_Type image_type);
      inline void set_scan_mode(Scan_Mode scan_mode);

      // --- Receiver Bandwidth --- //

      inline double  get_water_fat_shift(void) const;
      inline Time_us get_sampling_period(void) const;
      inline void    set_water_fat_shift(double water_fat_shift);
      inline void    set_sampling_period(Time_us sampling_period);

      // --- Convenience functions --- //
 
      const char *get_partial_fourier_name(Partial_Fourier_Type method) const;
      const char *get_orientation_name(Image_Orientation orientation) const;
      const char *get_image_type_name(Image_Type image_type) const;
      const char *get_scan_mode_name(Scan_Mode scan_mode) const;

      virtual void display_info(ostream& stream) const;

   protected:

      // --- VIO_Volume dimensions --- //

      double               _FOV[3];
      int                  _matrix_size[3];
      int                  _recon_size[3];
      double               _offset[3];
      double               _step[3];
      double               _slice_thickness;

      // --- Image Options --- //

      int                  _foldover_suppression;
      int                  _foldover_direction;
      Partial_Fourier_Type _partial_fourier_method;
      double               _scan_percentage;

      int                  _num_of_averages;
      Image_Orientation    _image_orientation;
      Image_Type           _image_type;
      Scan_Mode            _scan_mode;

      // --- Receiver Bandwidth --- //

      double               _water_fat_shift;

};

//--------------------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Pulse_Sequence::get_matrix_size
// Returns the length of the given dimension.
//--------------------------------------------------------------------------

inline 
int Pulse_Sequence::get_matrix_size(int i) const {

#ifdef DEBUG
   assert(i >= 0 && i < 3);
#endif

   return _matrix_size[i];
}

inline
void Pulse_Sequence::get_matrix_size(int matrix_size[]) const {
   matrix_size[0] = _matrix_size[0];
   matrix_size[1] = _matrix_size[1];
   matrix_size[2] = _matrix_size[2];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_reconstruction_size
//--------------------------------------------------------------------------

inline
int Pulse_Sequence::get_reconstruction_size(int i) const {
   return _recon_size[i];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_nslices
// Returns the number of slices for the pulse sequence.
//--------------------------------------------------------------------------

inline
int Pulse_Sequence::get_nslices(void) const {
   return _matrix_size[SLICE];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_nrows
// Returns the number of rows for the pulse sequence.
//--------------------------------------------------------------------------

inline
int Pulse_Sequence::get_nrows(void) const {
   return _matrix_size[ROW];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_ncols
// Returns the number of columns for the pulse sequence.
//--------------------------------------------------------------------------

inline
int Pulse_Sequence::get_ncols(void) const {
   return _matrix_size[COLUMN];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_fov
// Returns the field of view of the given dimension.
//--------------------------------------------------------------------------

inline 
double Pulse_Sequence::get_fov(int i) const {

#ifdef DEBUG
   assert(i >= 0 && i < 3);
#endif

   return _FOV[i];
}

inline
void Pulse_Sequence::get_fov(double fov[]) const {
   fov[0] = _FOV[0];
   fov[1] = _FOV[1];
   fov[2] = _FOV[2];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_voxel_offset
// Returns the offset of the given dimension.
//--------------------------------------------------------------------------

inline 
double Pulse_Sequence::get_voxel_offset(int i) const {

#ifdef DEBUG
   assert(i >= 0 && i < 3);
#endif

   return _offset[i];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_voxel_volume
// Returns the volume of a voxel (assume rectangular voxels).
//--------------------------------------------------------------------------

inline 
double Pulse_Sequence::get_voxel_volume(void) const {
   return get_slice_thickness() * get_voxel_step(ROW) * get_voxel_step(COLUMN);
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_slice_thickness
// Returns the pulse sequence slice thickness.
//--------------------------------------------------------------------------

inline
double Pulse_Sequence::get_slice_thickness(void) const {
   return _slice_thickness;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_fov
// Set the field of view for the pulse sequence.
//--------------------------------------------------------------------------

inline 
void Pulse_Sequence::set_fov(double fov[]) {

#ifdef DEBUG
   assert(fov != NULL);
#endif

   _FOV[0] = fov[0];
   _FOV[1] = fov[1];
   _FOV[2] = fov[2];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_matrix_size
// Set the matrix size for the pulse sequence.
//--------------------------------------------------------------------------

inline 
void Pulse_Sequence::set_matrix_size(int matrix_size[]) {

#ifdef DEBUG
   assert(matrix_size != NULL);
#endif

   _matrix_size[0] = matrix_size[0];
   _matrix_size[1] = matrix_size[1];
   _matrix_size[2] = matrix_size[2];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_voxel_offset
// Set the offset for the pulse sequence.
//--------------------------------------------------------------------------

inline 
void Pulse_Sequence::set_voxel_offset(double offset[]) {

#ifdef DEBUG
   assert(offset != NULL);
#endif

   _offset[0] = offset[0];
   _offset[1] = offset[1];
   _offset[2] = offset[2];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_slice_thickness
// Sets the slice thickness for the pulse sequence.
//--------------------------------------------------------------------------

inline
void Pulse_Sequence::set_slice_thickness(double slice_thickness) {
   _slice_thickness = slice_thickness;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_reconstruction_size
// Sets the size of the reconstruction matrix for the pulse sequence.
//--------------------------------------------------------------------------

inline
void Pulse_Sequence::set_reconstruction_size(int reconstruction_size[]) {
   _recon_size[0] = reconstruction_size[0];
   _recon_size[1] = reconstruction_size[1];
   _recon_size[2] = reconstruction_size[2];
}

//--------------------------------------------------------------------------
// Pulse_Sequence::uses_foldover_suppression
// Returns TRUE if foldover suppression is turned on.
//--------------------------------------------------------------------------

inline 
int Pulse_Sequence::uses_foldover_suppression(void) const {
   return _foldover_suppression;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::uses_partial_fourier
// Returns TRUE if a partial fourier method is turned on.
//--------------------------------------------------------------------------

inline 
int Pulse_Sequence::uses_partial_fourier(void) const {
   return (_partial_fourier_method != NO_PARTIAL_FOURIER);
}


//--------------------------------------------------------------------------
// Pulse_Sequence::get_partial_fourier_method
// Returns the type of partial fourier method being used.
//--------------------------------------------------------------------------

inline 
Partial_Fourier_Type Pulse_Sequence::get_partial_fourier_method(void) const {
   return _partial_fourier_method;
} 

//--------------------------------------------------------------------------
// Pulse_Sequence::get_scan_percentage
// Returns the scan percentage being used.
//--------------------------------------------------------------------------

inline 
double Pulse_Sequence::get_scan_percentage(void) const {
   return _scan_percentage;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_foldover_direction
// Returns the pulse sequence foldover direction.
//--------------------------------------------------------------------------

inline
int Pulse_Sequence::get_foldover_direction(void) const {
   return _foldover_direction;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_foldover_suppression
// Turns TRUE or FALSE foldover suppression.
//--------------------------------------------------------------------------

inline 
void Pulse_Sequence::set_foldover_suppression(int on) {

#ifdef DEBUG
   assert(on == TRUE || on == FALSE);
#endif

   _foldover_suppression = on;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_foldover_direction
// Sets the pulse sequence foldover direction.
//--------------------------------------------------------------------------

inline
void Pulse_Sequence::set_foldover_direction(int direction) {
   _foldover_direction = direction;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_partial_fourier_method
// Specifies a partial fourier method to use along with its scan percentage
// parameter.
//--------------------------------------------------------------------------

inline 
void Pulse_Sequence::set_partial_fourier_method(Partial_Fourier_Type method, 
                                                double scan_percentage) {
#ifdef DEBUG
   assert(scan_percentage >= 0.0);
   assert(scan_percentage <= 1.0);
#endif

   _partial_fourier_method = method;
   _scan_percentage        = scan_percentage;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_num_of_averages
// Returns the number of signal averages.
//--------------------------------------------------------------------------

inline 
int Pulse_Sequence::get_num_of_averages(void) const {
   return _num_of_averages;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_image_orientation
// Returns the image orientation.
//--------------------------------------------------------------------------

inline 
Image_Orientation Pulse_Sequence::get_image_orientation(void) const {
   return _image_orientation;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_image_type
// Returns the reconstructed image type.
//--------------------------------------------------------------------------

inline 
Image_Type Pulse_Sequence::get_image_type(void) const {
   return _image_type;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_scan_mode
// Returns the image scan mode.
//--------------------------------------------------------------------------

inline 
Scan_Mode Pulse_Sequence::get_scan_mode(void) const {
   return _scan_mode;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_num_of_averages
// Sets the number of signal averages.
//--------------------------------------------------------------------------

inline 
void Pulse_Sequence::set_num_of_averages(int num_of_averages) {

#ifdef DEBUG
   assert(num_of_averages >= 1);
#endif

   _num_of_averages = num_of_averages;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_image_orientation
// Sets the image orientation.
//--------------------------------------------------------------------------

inline 
void Pulse_Sequence::set_image_orientation(Image_Orientation orientation) {
   _image_orientation = orientation;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_image_type
// Sets the image type.
//--------------------------------------------------------------------------

inline 
void Pulse_Sequence::set_image_type(Image_Type image_type) {
   _image_type = image_type;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_scan_mode
// Sets the image scan mode.
//--------------------------------------------------------------------------

inline 
void Pulse_Sequence::set_scan_mode(Scan_Mode scan_mode) {
   _scan_mode = scan_mode;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_water_fat_shift
// Returns the water fat shift in pixels for the pulse sequence.
//--------------------------------------------------------------------------

inline 
double Pulse_Sequence::get_water_fat_shift(void) const {
   return _water_fat_shift;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::get_sampling_period
// Returns the sampling period in us for the pulse sequence.
//--------------------------------------------------------------------------

inline 
Time_us Pulse_Sequence::get_sampling_period(void) const {
   return (1.0E6 * _water_fat_shift) / 
          (WATER_FAT_SHIFT_HZ * get_matrix_size(_foldover_direction));
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_water_fat_shift
// Sets the water fat shift in pixels for the pulse sequence.
//--------------------------------------------------------------------------

inline 
void Pulse_Sequence::set_water_fat_shift(double water_fat_shift) {

#ifdef DEBUG
   assert(water_fat_shift > 0.0);
#endif

   _water_fat_shift = water_fat_shift;
}

//--------------------------------------------------------------------------
// Pulse_Sequence::set_sampling_period
// Sets the sampling period in us for the pulse sequence.
//--------------------------------------------------------------------------

void Pulse_Sequence::set_sampling_period(Time_us sampling_period) {

#ifdef DEBUG
   assert(sampling_period > 0.0);
   assert(get_matrix_size(READOUT) > 0);
#endif

   _water_fat_shift = sampling_period * 1.0E-6 *
                      WATER_FAT_SHIFT_HZ * get_matrix_size(READOUT);
}

#endif
