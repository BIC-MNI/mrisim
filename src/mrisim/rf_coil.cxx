//==========================================================================
// RF_COIL.CXX
// Base RF_coil class.
// Inherits from:  
// Base class to:  Image_SNR_Coil, Intrinsic_Coil, Percent_Coil
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/rf_coil.cxx,v 1.1 2003-05-30 16:43:11 bert Exp $
// $Log: rf_coil.cxx,v $
// Revision 1.1  2003-05-30 16:43:11  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  16:08:19  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:18:32  rkwan
// Release 2.5
//
// Revision 1.4  1996/01/17  18:02:51  rkwan
// Update for tx/rx coil modelling.
//
// Revision 1.3  1995/12/22  20:15:35  rkwan
// Update for percent coil and RF map features.
//
// Revision 1.2  1995/12/11  15:28:55  rkwan
// Doc update.
//
//==========================================================================

#include "rf_coil.h"

//---------------------------------------------------------------------------
// RF_Coil constructor
//---------------------------------------------------------------------------

RF_Coil::RF_Coil(long seed, const char *rx_map_file, const char *tx_map_file) {

   // --- Noise generation --- //

   // Initialize noise parameters
   _noise_mean     = 0.0;
   _noise_variance = 0.0;
   _coil_gain      = 1.0;

   // Set the random number generator random seed.
   // A random seed of 0 means to automatically generate a seed
   set_random_seed(seed);

   // --- Receive/Transmit Inhomogeneity maps --- //

   // Copy the receive inhomogeneity map file name
   if (rx_map_file == NULL){
      _rx_map_file = NULL;
   } else if ((_rx_map_file = strdup(rx_map_file)) == NULL){
      // Could not copy receive map file name
      cerr << "Out of memory." << endl;
      exit(EXIT_FAILURE);
   }

   // Copy the transmit inhomogeneity map file name
   if (tx_map_file == NULL){
      _tx_map_file = NULL;
   } else if ((_tx_map_file = strdup(tx_map_file)) == NULL){
      // Could not copy transmit map file
      cerr << "Out of memory." << endl;
      exit(EXIT_FAILURE);
   }

   // Open reception inhomogeneity map file
   if (this->uses_rx_map()) {
      _rx_map.open(_rx_map_file);
      _rx_map.set_default_float_icv();
      _rx_map.attach_icv();
   }

   // Open transmission inhomogeneity map file
   if (this->uses_tx_map()) {
      _tx_map.open(_tx_map_file);
      _tx_map.set_default_float_icv();
      _tx_map.attach_icv();
   }

}

//---------------------------------------------------------------------------
// RF_Coil Destructor
//---------------------------------------------------------------------------

RF_Coil::~RF_Coil() {

   // If a signal reception map was used, clean it up.
   if (this->uses_rx_map()){
      free(_rx_map_file);
      if (_rx_map.is_open()) {
         _rx_map.close();
      }
   }
   
   // If a signal transmission map was used, clean it up.
   if (this->uses_tx_map()){
      free(_tx_map_file);
      if (_tx_map.is_open()) {
         _tx_map.close();
      }
   }
}

//---------------------------------------------------------------------------
// RF_Coil::load_rx_map_slice
// Loads a slice of the receive coil map.
//---------------------------------------------------------------------------

void RF_Coil::load_rx_map_slice(int slice_num, Real_Slice& rx_slice) {

#ifdef DEBUG
   assert(this->uses_rx_map());
   assert(_rx_map.is_same_slice_size_as(rx_slice));
   assert(slice_num >= 0);
   assert(slice_num < _rx_map.get_nslices());
#endif

   _rx_map.load_slice(slice_num, (void *)rx_slice);
}

//---------------------------------------------------------------------------
// RF_Coil::load_tx_map_slice
// Loads a slice of the transmit coil map.
//---------------------------------------------------------------------------

void RF_Coil::load_tx_map_slice(int slice_num, Real_Slice& tx_slice) {

#ifdef DEBUG
   assert(this->uses_tx_map());
   assert(_tx_map.is_same_slice_size_as(tx_slice));
   assert(slice_num >= 0);
   assert(slice_num < _tx_map.get_nslices());
#endif

   _tx_map.load_slice(slice_num, (void *)tx_slice);
}

//---------------------------------------------------------------------------
// RF_Coil::get_rx_map_range
// Gets the maximum and minimum real ranges for a given slice of
// the receive coil volume.
//---------------------------------------------------------------------------

void RF_Coil::get_rx_map_range(int slice_num, double& min, double& max) const {
   long start[3] = {0, 0, 0};
   start[SLICE]  = slice_num;

   _rx_map.get_variable(MIimagemin, start, NC_DOUBLE, NULL, &min);
   _rx_map.get_variable(MIimagemax, start, NC_DOUBLE, NULL, &max);

}

//---------------------------------------------------------------------------
// RF_Coil::get_rx_map_range
// Gets the maximum and minimum real ranges for the entire receive coil
// volume.
//---------------------------------------------------------------------------

void RF_Coil::get_rx_map_range(double& min, double& max) const {

   double slice_min, slice_max;
   int    nslices = _rx_map.get_nslices();

   get_rx_map_range(--nslices, min, max);
 
   while (nslices-- > 0) {
      get_rx_map_range(nslices, slice_min, slice_max);
      if (slice_min < min) min = slice_min;
      if (slice_max > max) max = slice_max;
   }

}

//---------------------------------------------------------------------------
// RF_Coil::get_tx_map_range
// Gets the maximum and minimum real ranges for a given slice of
// the transmit coil volume.
//---------------------------------------------------------------------------

void RF_Coil::get_tx_map_range(int slice_num, double& min, double& max) const {
   long start[3] = {0, 0, 0};
   start[SLICE]  = slice_num;

   _tx_map.get_variable(MIimagemin, start, NC_DOUBLE, NULL, &min);
   _tx_map.get_variable(MIimagemax, start, NC_DOUBLE, NULL, &max);

}

//---------------------------------------------------------------------------
// RF_Coil::get_tx_map_range
// Gets the maximum and minimum real ranges for the entire transmit coil
// volume.
//---------------------------------------------------------------------------

void RF_Coil::get_tx_map_range(double& min, double& max) const {

   double slice_min, slice_max;
   int    nslices = _tx_map.get_nslices();

   get_tx_map_range(--nslices, min, max);
 
   while (nslices-- > 0) {
      get_rx_map_range(nslices, slice_min, slice_max);
      if (slice_min < min) min = slice_min;
      if (slice_max > max) max = slice_max;
   }

}

//---------------------------------------------------------------------------
// RF_Coil::display_coil_info
// Outputs coil information to an output stream.
//---------------------------------------------------------------------------

void RF_Coil::display_coil_info(ostream& stream) const {

   stream << "RF Coil Information:" << endl;
   stream << "--------------------" << endl << endl;
   stream << "Random seed:                " << _seed << endl;
   stream << "Noise mean:                 " << _noise_mean << endl;
   stream << "Noise variance:             " << _noise_variance << endl;
   stream << "Receive inhomogeneity map:  " 
          << (this->uses_rx_map() ? _rx_map_file : "NONE") << endl;
   stream << "Transmit inhomogeneity map: " 
          << (this->uses_tx_map() ? _tx_map_file : "NONE") << endl << endl;

}

//---------------------------------------------------------------------------
// RF_Coil::add_noise_to_real_image
// Returns a modified image with noise added according to the
// pulse sequence parameters.
//---------------------------------------------------------------------------

void RF_Coil::add_noise_to_real_image(Real_Slice& image_slice){
   double noise_real, noise_imag;

   if (_noise_variance == 0.0)
      return;  // do nothing and return

   // Add Gaussian noise to the image
   unsigned int m,n;
   for(m=0; m<image_slice.get_nrows(); m++){
      for(n=0; n<image_slice.get_ncols(); n++){
         _generate_gaussian_noise(noise_real, noise_imag);
         image_slice(m,n) = image_slice(m,n) + noise_real;
      } 
   }
}

//---------------------------------------------------------------------------
// RF_Coil::add_noise_to_modulus_image
// Returns a modified image with noise added according to the
// pulse sequence parameters.
//---------------------------------------------------------------------------

void RF_Coil::add_noise_to_modulus_image(Real_Slice& image_slice) {

   double noise_real, noise_imag;

   if (_noise_variance == 0.0)
      return;  // do nothing and return

   // Add Gaussian noise to inphase and quadrature channels
   unsigned int m,n;
   for(m=0; m<image_slice.get_nrows(); m++){
      for(n=0; n<image_slice.get_ncols(); n++){
         _generate_gaussian_noise(noise_real, noise_imag);
         image_slice(m,n) = hypotf( image_slice(m,n) + noise_real,
                                    noise_imag );
      } 
   }
   
}

//---------------------------------------------------------------------------
// RF_Coil::add_noise_to_complex_image
// Returns a modified complex image with noise added according to the
// pulse sequence parameters.
//---------------------------------------------------------------------------

void RF_Coil::add_noise_to_complex_image(Complex_Slice& image_slice) {

   double noise_real, noise_imag;

   if (_noise_variance == 0.0)
      return;  // do nothing and return

   // Add Gaussian noise to inphase and quadrature channels
   unsigned int m,n;
   for(m=0; m<image_slice.get_nrows(); m++){
      for(n=0; n<image_slice.get_ncols(); n++){
         _generate_gaussian_noise(noise_real, noise_imag);
         image_slice.real(m,n) += noise_real;
         image_slice.imag(m,n) += noise_imag;
      } 
   }

}

//---------------------------------------------------------------------------
// RF_Coil::add_noise_to_raw_slice
// Adds Gaussian noise to real and imaginary parts of a Complex_Matrix
// with variance according to the pulse sequence parameters.
//---------------------------------------------------------------------------

void RF_Coil::add_noise_to_raw_slice(Complex_Slice& raw_slice) {

   double noise_real, noise_imag;
   const double std_fft_scale = sqrt(raw_slice.get_nrows()*
                                     raw_slice.get_ncols());

   if (_noise_variance == 0.0)
      return;  // do nothing and return

   // Add Gaussian noise to inphase and quadrature channels
   unsigned int m,n;
   for(m=0; m<raw_slice.get_nrows(); m++){
      for(n=0; n<raw_slice.get_ncols(); n++){
         _generate_gaussian_noise(noise_real, noise_imag);
         raw_slice.real(m,n) += std_fft_scale*noise_real;
         raw_slice.imag(m,n) += std_fft_scale*noise_imag;
      } 
   }

}

//---------------------------------------------------------------------------
// RF_Coil::compute_variance
// Computes the noise variance for the current pulse sequence.
//---------------------------------------------------------------------------

void RF_Coil::_compute_variance(const Phantom &/*phantom*/){

   // Do nothing.
   // Noise variance is set by set_noise_variance() in base RF_Coil.

}

//---------------------------------------------------------------------------
// RF_Coil::_generate_gaussian_noise
// Returns two Gaussian distributed zero mean random numbers with
// the variance calculated for the current pulse sequence and scanner.
//
// Uses the Box-Muller method.
// Reference:  Press, William A. et al.  Numerical Recipes in C,
//             Second Edition, Cambridge Univ. Press, 1992, pp 288-290.
//---------------------------------------------------------------------------

void RF_Coil::_generate_gaussian_noise(double& n1, double& n2){
   double v1, v2;
   double R_squared, factor;

   do {

#ifdef NO_DRAND
      v1 = 2.0*rand()-1.0;
      v2 = 2.0*rand()-1.0;
#else
      v1 = 2.0*drand48()-1.0;
      v2 = 2.0*drand48()-1.0;
#endif
      R_squared = v1*v1 + v2*v2;
   } while (R_squared >= 1.0 || R_squared == 0.0);

   factor = sqrt(-2.0*log(R_squared)*_noise_variance/R_squared);
   n1 = v1*factor;
   n2 = v2*factor;
}
