//===========================================================================
// RF_TISSUE_PHANTOM.CXX
//
// R.Kwan
// April 29, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/rf_tissue_phantom.cxx,v 1.1 2003-05-30 16:43:12 bert Exp $
// $Log: rf_tissue_phantom.cxx,v $
// Revision 1.1  2003-05-30 16:43:12  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:19:24  rkwan
// Release 2.5
//
//===========================================================================

#include "rf_tissue_phantom.h"
#include <signal/quick_model.h>
#include <signal/isochromat_model.h>
#include <signal/fast_iso_model.h>

//---------------------------------------------------------------------------
// RF_Tissue_Phantom constructor
//---------------------------------------------------------------------------

RF_Tissue_Phantom::RF_Tissue_Phantom(unsigned int n_tissue_classes,
                                     unsigned int n_flip_angles) :
   Phantom(n_tissue_classes) {

#ifdef DEBUG
   assert(n_flip_angles > 0);
#endif

   // initialize table sizes
   // and allocate flip angle error tables

   _n_flip_angles   = n_flip_angles;
   _flip_error      = new double[_n_flip_angles];
   _flip_error_step = 0;
 
   // allocate signal intensity lookup tables
   // tables are indexed by tissue index

   _real_intensity = new double[_n_tissue_classes*_n_flip_angles];
   _imag_intensity = new double[_n_tissue_classes*_n_flip_angles];

   _no_error_real  = new double[_n_tissue_classes];
   _no_error_imag  = new double[_n_tissue_classes];

}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom destructor
//---------------------------------------------------------------------------

RF_Tissue_Phantom::~RF_Tissue_Phantom() {

   // Release allocated tables
   delete[] _flip_error;
   delete[] _real_intensity;
   delete[] _imag_intensity;
   delete[] _no_error_real;
   delete[] _no_error_imag;

}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::display_tissue_info
//---------------------------------------------------------------------------

void RF_Tissue_Phantom::display_tissue_info(ostream& stream) const {

   stream << "Tissue Simulation:" << endl; 
   stream << "------------------" << endl << endl;
   stream << "Number of tissues installed: " << _n_tissues_installed << endl;
   stream << "Number of flip angles used:  " << _n_flip_angles << endl
          << endl;
   stream << "ID:          Tissue Name     T1     T2    T2*     NH "
          << "Intensity (I,Q)" << endl
          << "--- -------------------- ------ ------ ------ ------ "
          << "------------------" 
          << endl;

   int tissue_index;
   for(tissue_index=0; tissue_index<_n_tissues_installed; tissue_index++){
      stream << setw(3) << (int)get_tissue_label(tissue_index) << " ";
      _tissue[tissue_index]->display_info(stream);
      stream << " " << setw(8) 
             << get_real_intensity(get_tissue_label(tissue_index))
             << " " << setw(8) 
             << get_imag_intensity(get_tissue_label(tissue_index))
             << endl;
   }
   stream << endl;

}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::apply_pulse_sequence
// Apply a pulse sequence to the tissue classes installed in the phantom.
//---------------------------------------------------------------------------

void RF_Tissue_Phantom::apply_pulse_sequence(Quick_Sequence *pseq) {
   unsigned int itissue, iflip;
   Quick_Model  *model;
   double       real, imag, mag;

      // Save pointer to current pulse sequence
   _pseq = pseq;

   // For each tissue class installed, apply the pulse sequence
   // and store the signal intensity.  Tissue classes which are
   // not used are assigned zero proton density so these are
   // assigned an intensity of 0.0

   for(itissue=0; itissue<_n_tissues_installed; itissue++){
      if (_tissue[itissue]->get_NH() != 0) {

         model = new Quick_Model(*_tissue[itissue]);

         pseq->apply(*model);
         real = _get_i_sample();
         imag = _get_q_sample();
         _no_error_real[itissue] = real;
         _no_error_imag[itissue] = imag;
         mag = hypot(real, imag);

         if (real >= _max_real) _max_real = real;
         if (real <  _min_real) _min_real = real;
         if (imag >= _max_imag) _max_imag = imag;
         if (imag <  _min_real) _min_imag = imag;
         if (mag >= _max_mag) _max_mag = mag;
         if (mag <  _min_mag) _min_mag = mag;

         if (uses_tx_map()){

            for (iflip=0; iflip<_n_flip_angles; iflip++){
               model->set_flip_error(_flip_error[iflip]);
               pseq->apply(*model);
               real = _get_i_sample();
               imag = _get_q_sample();
               _lookup_real_intensity(itissue, iflip) = real;
               _lookup_imag_intensity(itissue, iflip) = imag;
               mag  = hypot(real, imag);                                      

               if (real >= _max_real) _max_real = real;
               if (real <  _min_real) _min_real = real;
               if (imag >= _max_imag) _max_imag = imag;
               if (imag <  _min_real) _min_imag = imag;
               if (mag >= _max_mag) _max_mag = mag;
               if (mag <  _min_mag) _min_mag = mag;

            }
         }
 
         delete model;

      } else {

         _no_error_real[itissue] = 0.0;
         _no_error_imag[itissue] = 0.0;

         if (0.0 >= _max_real) _max_real = 0.0; 
         if (0.0 <  _min_real) _min_real = 0.0;
         if (0.0 >= _max_imag) _max_imag = 0.0;
         if (0.0 <  _min_real) _min_imag = 0.0; 
         if (0.0 >= _max_mag) _max_mag = 0.0;
         if (0.0 <  _min_mag) _min_mag = 0.0;

         if (uses_tx_map()){
            for (iflip=0; iflip<_n_flip_angles; iflip++){
               _lookup_real_intensity(itissue, iflip) = 0.0;
               _lookup_imag_intensity(itissue, iflip) = 0.0;

               if (0.0 >= _max_real) _max_real = 0.0; 
               if (0.0 <  _min_real) _min_real = 0.0;
               if (0.0 >= _max_imag) _max_imag = 0.0;
               if (0.0 <  _min_real) _min_imag = 0.0; 
               if (0.0 >= _max_mag) _max_mag = 0.0;
               if (0.0 <  _min_mag) _min_mag = 0.0;

            }
         }
      }
   }

}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::find_steady_state
// Apply a pulse sequence to the tissue classes installed in
// the phantom until a steady state is reached.
//---------------------------------------------------------------------------

void RF_Tissue_Phantom::find_steady_state(Custom_Sequence *pseq) {
   unsigned int          itissue, iflip;
   Fast_Isochromat_Model *model;
   double                real, imag, mag;

   // Save pointer to current pulse sequence
   _pseq = pseq;

   // For each tissue class installed, apply the pulse sequence
   // to steady state and store the signal intensity.  Tissue classes 
   // which are not used are assigned zero proton density so these are
   // assigned an intensity of 0.0
   
   for(itissue=0; itissue<_n_tissues_installed; itissue++){

      if (_tissue[itissue]->get_NH() != 0) {

         model = new Fast_Isochromat_Model(*_tissue[itissue]);

         pseq->initialize_sequence(*model);
         pseq->apply_to_steady_state(*model);
         real = _get_i_sample();
         imag = _get_q_sample();
         _no_error_real[itissue] = real;
         _no_error_imag[itissue] = imag;
         mag  = hypot(real, imag);

         if (real >= _max_real) _max_real = real;
         if (real <  _min_real) _min_real = real;
         if (imag >= _max_imag) _max_imag = imag;
         if (imag <  _min_real) _min_imag = imag;
         if (mag >= _max_mag) _max_mag = mag;
         if (mag <  _min_mag) _min_mag = mag;

         if (uses_tx_map()) {
            for (iflip=0; iflip<_n_flip_angles; iflip++){
               model->set_flip_error(_flip_error[iflip]);
               pseq->initialize_sequence(*model);
               pseq->apply_to_steady_state(*model);

               // Store the steady state magnetization
               real = _get_i_sample();
               imag = _get_q_sample();
               //_real_intensity[itissue*_n_flip_angles+iflip] = real;
               //_imag_intensity[itissue*_n_flip_angles+iflip] = imag;
               _lookup_real_intensity(itissue, iflip) = real;
               _lookup_imag_intensity(itissue, iflip) = imag;
               mag  = hypot(real, imag);

               if (real >= _max_real) _max_real = real;
               if (real <  _min_real) _min_real = real;
               if (imag >= _max_imag) _max_imag = imag;
               if (imag <  _min_real) _min_imag = imag;
               if (mag >= _max_mag) _max_mag = mag;
               if (mag <  _min_mag) _min_mag = mag;
            }
         }
         
         delete model;

      } else {

         _no_error_real[itissue] = 0.0;
         _no_error_imag[itissue] = 0.0;

         if (0.0 >= _max_real) _max_real = 0.0;
         if (0.0 <  _min_real) _min_real = 0.0;
         if (0.0 >= _max_imag) _max_imag = 0.0;
         if (0.0 <  _min_real) _min_imag = 0.0;
         if (0.0 >= _max_mag) _max_mag = 0.0;
         if (0.0 <  _min_mag) _min_mag = 0.0;

         if (uses_tx_map()) {

            for (iflip=0; iflip<_n_flip_angles; iflip++){
               //_real_intensity[itissue*_n_flip_angles+iflip] = 0.0;
               //_imag_intensity[itissue*_n_flip_angles+iflip] = 0.0;
               _lookup_real_intensity(itissue, iflip) = 0.0;
               _lookup_imag_intensity(itissue, iflip) = 0.0;

               if (0.0 >= _max_real) _max_real = 0.0;
               if (0.0 <  _min_real) _min_real = 0.0;
               if (0.0 >= _max_imag) _max_imag = 0.0;
               if (0.0 <  _min_real) _min_imag = 0.0;
               if (0.0 >= _max_mag) _max_mag = 0.0;
               if (0.0 <  _min_mag) _min_mag = 0.0;
            }

         }

      } 
   }

}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::_interp_real_intensity
// Returns the simulated intensity for Tissue tissue_label
//---------------------------------------------------------------------------

double RF_Tissue_Phantom::_interp_real_intensity(unsigned int index,
                                              double flip_error) const {

   // Interpolate flip angle error index
   unsigned int n = (int)floor((flip_error - _flip_error[0]) / 
                               _flip_error_step);

   const double X0 = _lookup_real_intensity(index, n);
   const double X1 = _lookup_real_intensity(index, n+1);
   return X0 + (X1 - X0) * (flip_error - _flip_error[n]) / _flip_error_step;
   
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::_interp_imag_intensity
// Returns the simulated intensity for Tissue tissue_label
//---------------------------------------------------------------------------

double RF_Tissue_Phantom::_interp_imag_intensity(unsigned int index,
                                              double flip_error) const {

   // Interpolate flip angle error index
   unsigned int n = (int)floor((flip_error - _flip_error[0]) /
                               _flip_error_step);

   const double X0 = _lookup_imag_intensity(index, n);
   const double X1 = _lookup_imag_intensity(index, n+1);
   return X0 + (X1 - X0) * (flip_error - _flip_error[n]) / _flip_error_step;
 
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::_interp_mag_intensity
// Returns the simulated intensity for Tissue tissue_label
//---------------------------------------------------------------------------

double RF_Tissue_Phantom::_interp_mag_intensity(unsigned int index,
                                             double flip_error) const {

   return hypot(_interp_real_intensity(index, flip_error),
                _interp_imag_intensity(index, flip_error));
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::get_signal_map
// Computes a nonuniform signal map for a tissue type given an RF coil
// transmit map.
//---------------------------------------------------------------------------

void RF_Tissue_Phantom::get_signal_map(unsigned int tissue_index,
                                       const Real_Slice& tx_map,
                                       Complex_Slice& signal_map) {
#ifdef DEBUG
   assert(tx_map.is_same_size_as(signal_map));
#endif

   unsigned int m, n;
   for (m=0; m<tx_map.get_nrows(); m++){
      for (n=0; n<tx_map.get_ncols(); n++){
         signal_map.real(m,n) = _interp_real_intensity(tissue_index, 
                                                       tx_map(m,n));
         signal_map.imag(m,n) = _interp_imag_intensity(tissue_index, 
                                                       tx_map(m,n));
      }
   }
}

//---------------------------------------------------------------------------
// RF_Tissue_Phantom::_setup_flip_errors
//---------------------------------------------------------------------------

void RF_Tissue_Phantom::_setup_flip_errors(double min, double max) {

   _flip_error_step = (max - min) / (_n_flip_angles - 1);
  
   int n;
   for(n=0; n<_n_flip_angles; n++){
      _flip_error[n] = n * _flip_error_step + min;
   } 

}

