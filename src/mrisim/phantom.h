#ifndef __PHANTOM_H
#define __PHANTOM_H

//===========================================================================
// PHANTOM.H
// Phantom class.
// Inherits from:
// Base class to:  Discrete_Phantom
//                 Fuzzy_Phantom
//                 RF_Phantom
//
// R.Kwan
// September 7, 1995
//
// (C) Copyright 1995, 1996 by R.Kwan
//===========================================================================

/*===========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/phantom.h,v 1.1 2003-05-30 16:43:11 bert Exp $
 * $Log: phantom.h,v $
 * Revision 1.1  2003-05-30 16:43:11  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  16:37:34  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:18:03  rkwan
 * Release 2.5
 *
 * Revision 1.3  1996/01/05  21:42:32  rkwan
 * Moved get_image_slice into Phantom, intensity calculations into MINCPhantom.
 *
 * Revision 1.3  1996/01/05  21:42:32  rkwan
 * Moved get_image_slice into Phantom, intensity calculations into MINCPhantom.
 *
 * Revision 1.2  1995/12/11  15:29:01  rkwan
 * Doc update.
 *
 *=========================================================================*/

#include <math.h>
#include <minc/mrimatrix.h>
#include <minc/mrilabel.h>
#include <minc/mriimage.h>
#include <minc/imincfile.h>
#include <minc/omincfile.h>
#include <minc/chirp.h>
#include <signal/quickseq.h>
#include <signal/customseq.h>
#include <signal/tissue.h>

typedef Label Tissue_Label;
typedef float Real_Scalar;
typedef MRI_Float_Matrix Real_Slice;
typedef MRI_FComplex_Matrix Complex_Slice;

#ifndef MAX_TISSUE_LABEL
#define MAX_TISSUE_LABEL MAX_LABEL
#endif

//---------------------------------------------------------------------------
// Phantom class
// Interface base class that stores spatial information about tissues and
// other characteristics in the phantom.   Implements interface to 
// generate simulated images.
//---------------------------------------------------------------------------

class Phantom {
   public:
      Phantom(unsigned int n_tissue_classes);

      virtual ~Phantom();

      // --- Tissue information --- //
      int install_tissue(Tissue* new_tissue, 
                         Tissue_Label tissue_label);
      inline unsigned int get_num_tissues(void) const;
      virtual void display_tissue_info(ostream& stream) const = 0;

      // --- Pulse Sequence simulation --- //
      virtual void apply_pulse_sequence(Quick_Sequence *pseq) = 0;
      virtual void find_steady_state(Custom_Sequence *pseq) = 0;

      // --- Image Formation --- //
      void ideal_lin_slice_select(double z_centre, 
                              double slice_thickness,
                              Real_Slice& output_slice);

      void ideal_nn_slice_select(double z_centre, 
                              double slice_thickness,
                              Real_Slice& output_slice);

      void initialize_chirp(unsigned int out_row_length,
                              unsigned int out_col_length,
                              double out_row_fov,
                              double out_col_fov);

      inline void generate_raw_data_slice(const Real_Slice& sim_slice,
                              Complex_Slice& raw_slice);

      void compute_partial_volume(const Real_Slice& sim_slice,
                              double row_step, double col_step,
                              double row_shift, double col_shift,
                              Real_Slice& pv_slice);  
   
      virtual void get_simulated_phantom_slice(int slice_num,
                              Complex_Slice& sim_slice) = 0;
      virtual void get_simulated_mag_phantom_slice(int slice_num,
                              Real_Slice& sim_slice) = 0;
      virtual void get_simulated_real_phantom_slice(int slice_num,
                              Real_Slice& sim_slice) = 0;

      inline void get_safe_simulated_phantom_slice(int slice_num,
                              Complex_Slice& sim_slice);
      inline void get_safe_simulated_mag_phantom_slice(int slice_num,
                              Real_Slice& sim_slice);
      inline void get_safe_simulated_real_phantom_slice(int slice_num,
                              Real_Slice& sim_slice);
                              
      // --- Access functions --- //
      virtual double get_real_intensity(Tissue_Label label) const = 0;
      virtual double get_imag_intensity(Tissue_Label label) const = 0;
      virtual double get_mag_intensity(Tissue_Label label) const = 0;

      inline Tissue_Label get_tissue_label(unsigned int index) const;
      inline unsigned int get_tissue_index(Tissue_Label tissue_label) const;
      inline int tissue_table_is_full(void) const;
      inline int tissue_label_is_valid(Tissue_Label tissue_label) const;

      inline void   get_max_intensity(double& real, double& imag) const;
      inline void   get_min_intensity(double& real, double& imag) const;
      inline double get_max_mag_intensity(void) const;
      inline double get_min_mag_intensity(void) const;
      inline double get_real_signal_gain(void) const;
      inline double get_imag_signal_gain(void) const;
      inline double get_mag_signal_gain(void) const;
  
      inline Pulse_Sequence &get_current_pulse_sequence(void) const;

      // --- VIO_Volume convenience functions --- //
      virtual int    is_same_slice_size_as(const MRI_Matrix& mat) const = 0;
      virtual int    get_nrows(void) const                              = 0;
      virtual int    get_ncols(void) const                              = 0;
      virtual int    get_nslices(void) const                            = 0;
      virtual void   get_volume_dimensions(int length[]) const          = 0;
      virtual void   get_volume_info(Volume_Info &volume_info) const    = 0;
      virtual double get_voxel_step(int n) const                        = 0;
      virtual double get_voxel_start(int n) const                       = 0;
      virtual void   display_volume_info(ostream& stream) const         = 0;
      virtual void   set_output_volume_info(O_MINC_File& output,
                                  const Volume_Info& vol_info, 
                                  const char *argstring = NULL) const   = 0;

      // --- Callback for Pulse Sequence interface --- //
      static void _save_sample(Vector_3D& v, void *obj);
   
   protected:

      // --- Internal member functions --- //
      void _direct_fourier_resample(const Complex_Slice& sim_slice,
                                    Complex_Slice& raw_slice);

      void _fft_fourier_resample(const Complex_Slice& sim_slice,
                                 Complex_Slice& raw_slice);

      void _chirp_fourier_resample(const Real_Slice& sim_slice,
                                   Complex_Slice& raw_slice);

      void _compute_chirp_weight(double weight[], 
                                 unsigned int in_length,
                                 unsigned int out_length,
                                 double in_fov, double out_fov);

      void _apply_weight(Real_Scalar in[], unsigned int stride, 
                         double weight[], unsigned int nelements);
 
      void _compensate_for_linear_kernel(Complex_Slice& raw_slice);

      // --- Pulse sequence simulation interface --- //
      inline double _get_i_sample(void) const;
      inline double _get_q_sample(void) const;

      double _min_real;   // minimum and maximum
      double _max_real;   // real channel signal
      double _min_imag;   // minimum and maximum
      double _max_imag;   // imaginary channel signal
      double _min_mag;    // minimum and maximum
      double _max_mag;    // magnitude signal

      Pulse_Sequence *_pseq;          // Current pulse sequence used to
                                      // compute signal intensities

      // --- Internal data structures --- //
      unsigned int  _n_tissue_classes;
      unsigned int  _n_tissues_installed;

      unsigned int  *_tissue_index;   // label -> index translation table
      Tissue_Label  *_tissue_label;   // index -> label translation table
      Tissue        **_tissue;        // Tissue tables (stored by index)

   private:

      // --- Pulse sequence simulation interface --- //
      static double       _I_sample;   // In-phase and
      static double       _Q_sample;   // Quadrature channels

      // --- Fourier Resampling --- //
      Chirp_Algorithm     *row_chirp;
      Chirp_Algorithm     *col_chirp;
      double              *row_weight;
      double              *col_weight;
      Complex_Slice       *tmp_slice;

};

//---------------------------------------------------------------------------
// _sinc
// Computes the sinc function sin(PI x) / (PI x) handling the singularity at
// (x == 0.0) correctly.
//---------------------------------------------------------------------------

inline
double _sinc(double x) {
   return ((x!=0.0) ? sin(M_PI*x)/(M_PI*x) : 1.0);
}

//---------------------------------------------------------------------------
// Inline member functions
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Phantom::get_num_tissues
// Returns the number of tissues installed in the phantom.
//---------------------------------------------------------------------------

inline
unsigned int Phantom::get_num_tissues(void) const {
   return _n_tissues_installed;
}

//---------------------------------------------------------------------------
// Phantom::generate_raw_data_slice
// Generates a raw data slice by resampling the phantom Fourier data.
//---------------------------------------------------------------------------

inline
void Phantom::generate_raw_data_slice(const Real_Slice& sim_slice,
                                      Complex_Slice& raw_slice) {

   _chirp_fourier_resample(sim_slice, raw_slice);
   //_compensate_for_linear_kernel(raw_slice);
}

//---------------------------------------------------------------------------
// Phantom::get_safe_simulated_phantom_slice
// Generates a coloured phantom slice from labelled phantom data.
// Checks to see if the slice number is valid, if not returns a
// simulated slice of zeros.
//---------------------------------------------------------------------------

inline
void Phantom::get_safe_simulated_phantom_slice(int slice_num,
                                            Complex_Slice& sim_slice) {

   if ((slice_num >= 0) && (slice_num < get_nslices())) {
      get_simulated_phantom_slice(slice_num, sim_slice);
   } else {
      sim_slice.zeros();
   }

}

//---------------------------------------------------------------------------
// Phantom::get_safe_simulated_mag_phantom_slice
// Generates a coloured phantom slice from labelled phantom data.
// Checks to see if the slice number is valid, if not returns a
// simulated slice of zeros.
//---------------------------------------------------------------------------

inline
void Phantom::get_safe_simulated_mag_phantom_slice(int slice_num,
                                            Real_Slice& sim_slice) {

   if ((slice_num >= 0) && (slice_num < get_nslices())) {
      get_simulated_mag_phantom_slice(slice_num, sim_slice);
   } else {
      sim_slice.zeros();
   }

}

//---------------------------------------------------------------------------
// Phantom::get_safe_simulated_real_phantom_slice
// Generates a coloured phantom slice from labelled phantom data.
// Checks to see if the slice number is valid, if not returns a
// simulated slice of zeros.
//---------------------------------------------------------------------------

inline
void Phantom::get_safe_simulated_real_phantom_slice(int slice_num,
                                            Real_Slice& sim_slice) {

   if ((slice_num >= 0) && (slice_num < get_nslices())) {
      get_simulated_real_phantom_slice(slice_num, sim_slice);
   } else {
      sim_slice.zeros();
   }

}

//---------------------------------------------------------------------------
// Phantom::get_tissue_label
// Lookup table index to Tissue_Label translation.
//---------------------------------------------------------------------------

inline
Tissue_Label Phantom::get_tissue_label(unsigned int index) const {
   return _tissue_label[index];
}

//---------------------------------------------------------------------------
// Phantom::get_tissue_index
// Tissue_Label to lookup table index translation.
//---------------------------------------------------------------------------

inline
unsigned int Phantom::get_tissue_index(Tissue_Label tissue_label) const{
  return _tissue_index[tissue_label];
}

//---------------------------------------------------------------------------
// Phantom::tissue_table_is_full
// Returns TRUE if the internal tissue table is full.
//---------------------------------------------------------------------------

inline
int Phantom::tissue_table_is_full(void) const {
   return (_n_tissues_installed == _n_tissue_classes);
}

//---------------------------------------------------------------------------
// Phantom::tissue_label_is_valid
// Returns TRUE if the tissue_label has been registered with the Phantom.
//---------------------------------------------------------------------------

inline
int Phantom::tissue_label_is_valid(Tissue_Label tissue_label) const {
   return ((tissue_label == 0) || (_tissue_index[tissue_label] != 0));
}

//---------------------------------------------------------------------------
// Tissue_Phantom::get_max_intensity
// Returns the maximum real and imaginary signal intensities.
//---------------------------------------------------------------------------

inline 
void Phantom::get_max_intensity(double& real, double& imag) const {
   real = _max_real;
   imag = _max_imag;
}

//---------------------------------------------------------------------------
// Tissue_Phantom::get_min_intensity
// Returns the minimum real and imaginary signal intensities.
//---------------------------------------------------------------------------

inline 
void Phantom::get_min_intensity(double& real, double& imag) const {
   real = _min_real;
   imag = _min_imag;
}

//---------------------------------------------------------------------------
// Tissue_Phantom::get_max_mag_intensity
// Returns the maximum magnitude signal intensity.
//---------------------------------------------------------------------------

inline 
double Phantom::get_max_mag_intensity(void) const {
   return _max_mag;
}

//---------------------------------------------------------------------------
// Tissue_Phantom::get_min_mag_intensity
// Returns the minimum magnitude signal intensity.
//---------------------------------------------------------------------------

inline 
double Phantom::get_min_mag_intensity(void) const {
   return _min_mag;
}

//---------------------------------------------------------------------------
// Phantom::get_real_signal_gain
// Returns the signal gain for the real channel.
//---------------------------------------------------------------------------

inline  
double Phantom::get_real_signal_gain(void) const {
   double gain;
   if (fabs(_max_real) > fabs(_min_real)){
      gain = 4095.0/fabs(_max_real);
   } else {
      gain = 4095.0/fabs(_min_real);
   }
   return gain;
}

//---------------------------------------------------------------------------
// Phantom::get_imag_signal_gain
// Returns the signal gain for the imaginary channel.
//---------------------------------------------------------------------------

inline
double Phantom::get_imag_signal_gain(void) const {
   double gain;
   if (fabs(_max_imag) > fabs(_min_imag)){
      gain = 4095.0/fabs(_max_imag);
   } else {
      gain = 4095.0/fabs(_min_imag);
   }
   return gain;
}

//---------------------------------------------------------------------------
// Phantom::get_mag_signal_gain
// Returns the signal gain for the magnitude of real and imaginary channels.
//---------------------------------------------------------------------------

inline
double Phantom::get_mag_signal_gain(void) const {
   double gain;
   if (fabs(_max_mag) > fabs(_min_mag)){
      gain = 4095.0/fabs(_max_mag);
   } else {
      gain = 4095.0/fabs(_min_mag);
   }
   return gain;
}

//---------------------------------------------------------------------------
// Phantom::get_current_pulse_sequence
// Returns a reference to the current pulse sequence.
//---------------------------------------------------------------------------

inline 
Pulse_Sequence &Phantom::get_current_pulse_sequence(void) const {
   return *_pseq;
}

//---------------------------------------------------------------------------
// Phantom::_save_sample
// Callback function to register with the Sample object.
// Saves the sample value in the protected Phantom member _sample.
//---------------------------------------------------------------------------

inline
void Phantom::_save_sample(Vector_3D& v, void *){
   Phantom::_I_sample = v[Y_AXIS];
   Phantom::_Q_sample = v[X_AXIS];
}

//---------------------------------------------------------------------------
// Phantom::_get_i_sample
// Return the computed in-phase sample.
//---------------------------------------------------------------------------

inline
double Phantom::_get_i_sample(void) const {
   return Phantom::_I_sample;
}

//---------------------------------------------------------------------------
// Phantom::_get_q_sample
// Return the computed quadrature sample.
//---------------------------------------------------------------------------

inline
double Phantom::_get_q_sample(void) const {
   return Phantom::_Q_sample;
}

#endif

