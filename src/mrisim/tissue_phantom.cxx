//===========================================================================
// TISSUE_PHANTOM.CXX
//
// R.Kwan
// April 29, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/tissue_phantom.cxx,v 1.1 2003-05-30 16:43:12 bert Exp $
// $Log: tissue_phantom.cxx,v $
// Revision 1.1  2003-05-30 16:43:12  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:19:48  rkwan
// Release 2.5
//
//===========================================================================

#include "tissue_phantom.h"
#include <signal/quick_model.h>
#include <signal/isochromat_model.h>
#include <signal/fast_iso_model.h>

//---------------------------------------------------------------------------
// Tissue_Phantom constructor
//---------------------------------------------------------------------------

Tissue_Phantom::Tissue_Phantom(unsigned int n_tissue_classes) :
   Phantom(n_tissue_classes) {

#ifdef DEBUG
   assert(n_tissue_classes > 0);
#endif

   // allocate signal intensity lookup tables
   // tables are indexed by tissue index

   _real_intensity = new double[_n_tissue_classes];
   _imag_intensity = new double[_n_tissue_classes];

}

//---------------------------------------------------------------------------
// Tissue_Phantom destructor
//---------------------------------------------------------------------------

Tissue_Phantom::~Tissue_Phantom() {

   // Release allocated tables
   delete[] _real_intensity;
   delete[] _imag_intensity;

}

//---------------------------------------------------------------------------
// Tissue_Phantom::display_tissue_info
//---------------------------------------------------------------------------

void Tissue_Phantom::display_tissue_info(ostream& stream) const {

   stream << "Tissue Simulation:" << endl; 
   stream << "------------------" << endl << endl;
   stream << "Number of tissues installed: " << _n_tissues_installed 
          << endl << endl;
   stream << "ID:          Tissue Name     T1     T2    T2*     NH "
          << "Intensity (I,Q)" << endl
          << "--- -------------------- ------ ------ ------ ------ "
          << "------------------" 
          << endl;

   int tissue_index;
   for(tissue_index=0; tissue_index<_n_tissues_installed; tissue_index++){
      stream << setw(3) << (int)get_tissue_label(tissue_index) << " ";
      _tissue[tissue_index]->display_info(stream);
      stream << " " << setw(8) << _lookup_real_intensity(tissue_index) 
             << " " << setw(8) << _lookup_imag_intensity(tissue_index)
             << endl;
   }
   stream << endl;

}

//---------------------------------------------------------------------------
// Tissue_Phantom::apply_pulse_sequence
// Apply a pulse sequence to the tissue classes installed in the phantom.
//---------------------------------------------------------------------------

void Tissue_Phantom::apply_pulse_sequence(Quick_Sequence *pseq) {
   unsigned int itissue;
   Quick_Model  *model;
   double       real, imag, mag;

      // Save pointer to current pulse sequence
   _pseq = pseq;

   // For each tissue class installed, apply the pulse sequence
   // and store the signal intensity.  Tissue classes which are
   // not used are assigned zero proton density so these are
   // assigned an intensity of 0.0

   for(itissue=0; itissue<_n_tissues_installed; itissue++){

      if (_tissue[itissue]->get_NH() != 0){

         model = new Quick_Model(*_tissue[itissue]);
         pseq->apply(*model);

         real = _real_intensity[itissue] = _get_i_sample();
         imag = _imag_intensity[itissue] = _get_q_sample();
         mag  = hypot(real, imag);

         delete model;

      } else {

         real = _real_intensity[itissue] = 0.0;
         imag = _imag_intensity[itissue] = 0.0;
         mag  = 0.0;

      }

      // Save min / max intensity
      if (real >= _max_real) _max_real = real;
      if (real <  _min_real) _min_real = real;
      if (imag >= _max_imag) _max_imag = imag;
      if (imag <  _min_real) _min_imag = imag;
      if (mag >= _max_mag) _max_mag = mag;
      if (mag <  _min_mag) _min_mag = mag;
   }

}

//---------------------------------------------------------------------------
// Tissue_Phantom::find_steady_state
// Apply a pulse sequence to the tissue classes installed in
// the phantom until a steady state is reached.
//---------------------------------------------------------------------------

void Tissue_Phantom::find_steady_state(Custom_Sequence *pseq) {
   int                   itissue;
   Fast_Isochromat_Model *model;
   double                real, imag, mag;

   // Save pointer to current pulse sequence
   _pseq = pseq;

   // For each tissue class installed, apply the pulse sequence
   // to steady state and store the signal intensity.  Tissue classes 
   // which are not used are assigned zero proton density so these are
   // assigned an intensity of 0.0
   
   for(itissue=0; itissue<_n_tissues_installed; itissue++){
      if (_tissue[itissue]->get_NH() != 0){

         // Instantiate a new magnetization system and run
         // the first repetition of the pulse sequence.
         model = new Fast_Isochromat_Model(*_tissue[itissue]);
         pseq->initialize_sequence(*model);
         pseq->apply_to_steady_state(*model);

         // Store the steady state magnetization
         real = _real_intensity[itissue] = _get_i_sample();
         imag = _imag_intensity[itissue] = _get_q_sample();
         mag  = hypot(real, imag);
         delete model;

      } else {

         real = _real_intensity[itissue] = 0.0;
         imag = _imag_intensity[itissue] = 0.0;
         mag  = 0.0;

      } 

      // Save min / max intensity
      if (real >= _max_real) _max_real = real;
      if (real <  _min_real) _min_real = real;
      if (imag >= _max_imag) _max_imag = imag;
      if (imag <  _min_real) _min_imag = imag;
      if (mag >= _max_mag) _max_mag = mag;
      if (mag <  _min_mag) _min_mag = mag;

   }

}
