#ifndef __RF_COIL_H
#define __RF_COIL_H

//==========================================================================
// RF_COIL.H
// Base RF_coil class.
// Inherits from:  
// Base class to:  Image_SNR_Coil, Intrinsic_Coil, Percent_Coil
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/rf_coil.h,v 1.1 2003-05-30 16:43:12 bert Exp $
 * $Log: rf_coil.h,v $
 * Revision 1.1  2003-05-30 16:43:12  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  16:08:10  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:18:27  rkwan
 * Release 2.5
 *
 * Revision 1.4  1996/01/17  18:02:57  rkwan
 * Update for tx/rx coil modelling.
 *
 * Revision 1.3  1995/12/22  20:15:30  rkwan
 * Update for percent coil and RF map features.
 *
 * Revision 1.2  1995/12/11  15:28:48  rkwan
 * Doc update.
 *
 *========================================================================*/

#include "phantom.h"
#include <minc/mrimatrix.h>
#include <minc/mriimage.h>
#include <minc/imincfile.h>

// For local time based random seeds
#include <sys/types.h>
#include <time.h>

//--------------------------------------------------------------------------
// RF_Coil class
//--------------------------------------------------------------------------

class RF_Coil {
   public:
      RF_Coil(long seed = 0, 
              const char *rx_map_file = NULL,
              const char *tx_map_file = NULL);

      virtual ~RF_Coil();

      // --- Random number generation --- //

      inline long set_random_seed(long seed = 0);
      inline long get_random_seed(void) const;

      inline void set_noise_variance(double variance);
      inline void set_noise_mean(double mean);
      inline double get_noise_variance(void) const;
      inline double get_noise_mean(void) const;

      // --- Coil inhomogeneity --- //

      inline int uses_rx_map(void) const;
      inline int uses_tx_map(void) const;

      inline int rx_map_is_same_size_as(const MINC_File&) const;
      inline int tx_map_is_same_size_as(const MINC_File&) const;
      inline int rx_map_has_same_dimension_info_as(const MINC_File&) const;
      inline int tx_map_has_same_dimension_info_as(const MINC_File&) const;
      inline int rx_map_has_same_dimension_names_as(const MINC_File&) const;
      inline int tx_map_has_same_dimension_names_as(const MINC_File&) const;

      void load_rx_map_slice(int slice_num, Real_Slice& rx_slice);
      void load_tx_map_slice(int slice_num, Real_Slice& tx_slice);

      void get_rx_map_range(int slice_num, double& min, double& max) const;
      void get_tx_map_range(int slice_num, double& min, double& max) const;
      void get_rx_map_range(double& min, double& max) const;
      void get_tx_map_range(double& min, double& max) const;
     
      // --- Convenience functions --- //

      virtual void display_coil_info(ostream& stream) const;

      // --- Noise simulation --- //

      void add_noise_to_real_image(Real_Slice& image_slice);
      void add_noise_to_modulus_image(Real_Slice& image_slice);

      void add_noise_to_complex_image(Complex_Slice& image_slice);
      void add_noise_to_raw_slice(Complex_Slice& raw_slice);

      virtual void _compute_variance(const Phantom &phantom);
      
   protected:

      // --- Internal data structures --- //

      double _noise_mean;       // mean of the noise
      double _noise_variance;   // variance of the noise
      long   _seed;             // random number generator seed
      double _coil_gain;        // real gain level for the coil

   private:

      // --- Internal member functions --- //
      void _generate_gaussian_noise(double& n1, double& n2);

      // --- Signal inhomogeneity maps --- //

      char   *_rx_map_file;     // Inhomogeneity map for signal
      I_MINC_File _rx_map;      // reception

      char   *_tx_map_file;     // Inhomogeneity map for signal
      I_MINC_File _tx_map;      // transmission

};

//--------------------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------------------

//---------------------------------------------------------------------------
// RF_Coil::set_random_seed
// Sets the seed for the random number generator and returns the seed
// value.   If a seed of 0 is given, sets the seed using the current
// local time.
//---------------------------------------------------------------------------

inline
long RF_Coil::set_random_seed(long seed){

   // If a seed of 0 is used, a seed will be selected based on
   // the current local time.
   if (seed == 0) {
      (void)time((time_t *)&seed);
   }
   _seed = seed;

#ifdef NO_DRAND
   srand(_seed);
#else
   srand48(_seed);
#endif

   return _seed;
}

//---------------------------------------------------------------------------
// RF_Coil::get_random_seed(void)
// Gets the seed for the random number generator.
//---------------------------------------------------------------------------

inline
long RF_Coil::get_random_seed(void) const {
   return _seed;
}

//---------------------------------------------------------------------------
// RF_Coil::set_noise_variance
// Access function to set the noise variance of the RF_Coil
//---------------------------------------------------------------------------

inline
void RF_Coil::set_noise_variance(double variance){

#ifdef DEBUG
   // Noise variance must be non negative.
   assert(_noise_variance >= 0.0);
#endif

   _noise_variance = variance;
}

//---------------------------------------------------------------------------
// RF_Coil::set_noise_mean
// Access function to set the noise mean of the RF_Coil
//---------------------------------------------------------------------------

inline
void RF_Coil::set_noise_mean(double mean){
   _noise_mean = mean;
}

//---------------------------------------------------------------------------
// RF_Coil::get_noise_variance
// Access function to get the noise variance of the RF_Coil
//---------------------------------------------------------------------------

inline
double RF_Coil::get_noise_variance(void) const {
   return _noise_variance;
}

//---------------------------------------------------------------------------
// RF_Coil::get_noise_mean
// Access function to get the noise mean of the RF_Coil
//---------------------------------------------------------------------------

inline
double RF_Coil::get_noise_mean(void) const {
   return _noise_mean;
}

//--------------------------------------------------------------------------
// RF_Coil::uses_rx_map
// Returns TRUE if the RF_Coil uses an inhomogeneity map for reception.
//--------------------------------------------------------------------------

inline
int RF_Coil::uses_rx_map(void) const {
   return (_rx_map_file != NULL);
}

//--------------------------------------------------------------------------
// RF_Coil::uses_tx_map
// Returns TRUE if the RF_Coil uses an inhomogeneity map for transmission.
//--------------------------------------------------------------------------

inline
int RF_Coil::uses_tx_map(void) const {
   return (_tx_map_file != NULL);
}

//--------------------------------------------------------------------------
// RF_Coil::rx_map_is_same_size_as
// Returns TRUE if minc_file is the same size as the receive map.
//--------------------------------------------------------------------------

inline
int RF_Coil::rx_map_is_same_size_as(const MINC_File &minc_file) const {
   return _rx_map.is_same_size_as(minc_file);
}

//--------------------------------------------------------------------------
// RF_Coil::tx_map_is_same_size_as
// Returns TRUE if minc_file is the same size as the transmit map.
//--------------------------------------------------------------------------

inline
int RF_Coil::tx_map_is_same_size_as(const MINC_File &minc_file) const {
   return _tx_map.is_same_size_as(minc_file);
}

//---------------------------------------------------------------------------
// RF_Coil::rx_map_has_same_dimension_info_as
// Returns TRUE if minc_file has the same start/step as the receive map.
//---------------------------------------------------------------------------

inline
int RF_Coil::rx_map_has_same_dimension_info_as(
      const MINC_File &minc_file) const {
   return _rx_map.has_same_dimension_info_as(minc_file);
}

//---------------------------------------------------------------------------
// RF_Coil::tx_map_has_same_dimension_info_as
// Returns TRUE if minc_file has the same start/step as the transmit map.
//---------------------------------------------------------------------------

inline
int RF_Coil::tx_map_has_same_dimension_info_as(
      const MINC_File &minc_file) const {
   return _tx_map.has_same_dimension_info_as(minc_file);
}

//---------------------------------------------------------------------------
// RF_Coil::rx_map_has_same_dimension_names_as
// Returns TRUE if minc_file has the same dimension names as the 
// receive map.
//---------------------------------------------------------------------------

inline
int RF_Coil::rx_map_has_same_dimension_names_as(
      const MINC_File &minc_file) const {
   return _rx_map.has_same_dimension_names_as(minc_file);
}

//---------------------------------------------------------------------------
// RF_Coil::tx_map_has_same_dimension_names_as
// Returns TRUE if minc_file has the same dimension names as the 
// transmit map.
//---------------------------------------------------------------------------

inline
int RF_Coil::tx_map_has_same_dimension_names_as(
      const MINC_File &minc_file) const {
   return _tx_map.has_same_dimension_names_as(minc_file);
}

#endif
