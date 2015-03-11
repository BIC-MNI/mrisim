//===========================================================================
// PHANTOM.CXX 
//
// R.Kwan
// September 7, 1995
//
// (C) Copyright 1995 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/phantom.cxx,v 1.1 2003-05-30 16:43:11 bert Exp $
// $Log: phantom.cxx,v $
// Revision 1.1  2003-05-30 16:43:11  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  16:37:53  rkwan
// Release 3.1 update.
//
// Revision 2.6  1996/05/29  19:14:44  rkwan
// GNU warning fix.
//
// Revision 2.5  1996/05/29  16:18:08  rkwan
// Release 2.5
//
// Revision 1.5  1996/01/17  18:02:41  rkwan
// Update for tx/rx coil modelling.
//
// Revision 1.4  1996/01/05  21:43:17  rkwan
// Moved get_image_slice into Phantom, intensity calculations into MINCPhantom.
//
// Revision 1.3  1995/12/22  20:16:43  rkwan
// Update for percent coil and RF map features.
//
// Revision 1.2  1995/12/11  15:29:06  rkwan
// Doc update.
//
//===========================================================================

#include "phantom.h"
#include <assert.h>
#include <float.h>

#include <signal/quick_model.h>
#include <signal/isochromat_model.h>
#include <signal/fast_iso_model.h>

#ifndef SLICE
#define SLICE 0
#define ROW   1
#define COLUMN 2
#endif

//---------------------------------------------------------------------------
// Static data members
//---------------------------------------------------------------------------

double Phantom::_I_sample;
double Phantom::_Q_sample;

//---------------------------------------------------------------------------
// Phantom constructor
//---------------------------------------------------------------------------

Phantom::Phantom(unsigned int n_tissue_classes) {

#ifdef DEBUG
   assert(n_tissue_classes > 0);
#endif

   // --- Allocate tissue lookup tables --- //

   // initialize table sizes
   // the first table entry is reserved for background by setting
   // _n_tissues_installed to one.

   _n_tissue_classes    = n_tissue_classes;
   _n_tissues_installed = 1;

   // index <-> label translation tables:
   // _tissue_index provides tissue label -> tissue index translation
   // _tissue_label provides tissue index -> tissue label translation

   _tissue_index = new unsigned int[MAX_TISSUE_LABEL];
   _tissue_label = new Tissue_Label[_n_tissue_classes];

   // allocate tissue lookup tables
   // tables are indexed by tissue index

   _tissue       = new Tissue *[_n_tissue_classes];

   // --- Initialize tissue lookup tables --- //

   int n;
   for(n=0; n<MAX_TISSUE_LABEL; n++){

      // Any undefined labels will point to tissue_label == 0
      // which is reserved for 'background'
      _tissue_index[n] = 0;

   }
   for(n=0; n<_n_tissue_classes; n++){
      _tissue[n] = (Tissue *)NULL;
      _tissue_label[n] = 0;
   }

   // Initialize signal max/min
   _max_real = DBL_MIN;
   _max_imag = DBL_MIN;
   _max_mag  = DBL_MIN;
   _min_real = DBL_MAX;
   _min_imag = DBL_MAX;
   _min_mag  = DBL_MAX;

   // Fourier resampling
   row_chirp  = NULL;
   col_chirp  = NULL;
   row_weight = NULL;
   col_weight = NULL;
   tmp_slice  = NULL;
}

//---------------------------------------------------------------------------
// Phantom destructor
//---------------------------------------------------------------------------

Phantom::~Phantom(){

   delete[] _tissue_index;
   delete[] _tissue_label;

   // Tissues installed in the phantom become the property of the
   // phantom and must be deleted when the phantom is deleted
   int n;
   for(n=0; n<_n_tissues_installed; n++){
      delete _tissue[n];
   }
   delete[] _tissue;

   // Clean up Fourier resampling temporaries

   if (row_chirp != NULL) delete row_chirp;
   if (col_chirp != NULL) delete col_chirp;
   if (row_weight != NULL) delete [] row_weight;
   if (col_weight != NULL) delete [] col_weight;
   if (tmp_slice != NULL) delete tmp_slice;

}

//---------------------------------------------------------------------------
// Phantom::install_tissue
//
// Installs tissue classes in the phantom and associates them with
// tissue labels in the labelled volume.
// Tissues are specified by a pointer to a tissue object.  These objects
// become the property of the Phantom and are destroyed with the phantom.
// Returns TRUE if the tissue is installed correctly, and FALSE
// if there is no room for the tissue to be installed.
//---------------------------------------------------------------------------

int Phantom::install_tissue(Tissue* new_tissue,
                            Tissue_Label tissue_label) {

   unsigned int end = _n_tissues_installed;
   int full = FALSE;

   // Check if the tissue table is full

   if (tissue_table_is_full()){

      full = TRUE;

   } else if (tissue_label == 0) {

      _tissue_index[0] = 0;
      _tissue_label[0] = 0;
      _tissue[0]       = new_tissue;

   } else {

      // Add new tissue to end of the tissue table
      // and update lookup table indices

      _tissue_index[tissue_label] = end;
      _tissue_label[end]          = tissue_label;
      _tissue[end]                = new_tissue;

      _n_tissues_installed++;

   }

   return full;

}

//---------------------------------------------------------------------------
// Phantom::ideal_lin_slice_select
// 
// Performs an ideal rectangular slice selection, for a given slice_thickness
// with the slice centred at z_centre.  Stores the result in output_slice.
// Phantom data is assumed to be linearly interpolated between samples
// and the slice selection profile is rectangular.
//
// Slice selection weighting is performed by integrating over the portion
// of the triangular interpolation kernel that falls within each phantom
// slice boundary.   
//---------------------------------------------------------------------------

void Phantom::ideal_lin_slice_select(double z_centre,
                                     double slice_thickness,
                                     Real_Slice& output_slice){

#ifdef DEBUG
   // Check that output_slice is the right size
   assert(this->is_same_slice_size_as(output_slice));
   assert(slice_thickness > 0.0);
#endif

   // Get the slice separation thickness in the phantom.
   double phantom_thickness = get_voxel_step(SLICE);

#ifdef DEBUG
   assert(phantom_thickness > 0.0);
#endif

   // Compute slice coordinate limits
   // Phantom z-coordinate starts at 0.
   double z_min = (z_centre - 0.5*slice_thickness)/phantom_thickness;
   double z_max = (z_centre + 0.5*slice_thickness)/phantom_thickness;

   int z_start = (int)floor(z_min);
   int z_stop  = (int)floor(z_max);

#ifdef TRACE
   cout << endl;
   cout << "z_centre: " << z_centre << " slice thickness: "
        << slice_thickness << endl;
   cout << "z_min: " << z_min << " z_max: " << z_max << endl;
   cout << "z_start: " << z_start << " z_stop: " << z_stop << endl;
#endif

   // Buffer two phantom slices at a time to compute a single output slice.  
   // Clear the output slice and create two buffer slices of the same size.
   // Set up a buffer of pointers to the slices and a temporary pointer
   // for swapping the buffers.

   output_slice.zeros();
   Real_Slice *buffer[2];
   buffer[0] = new Real_Slice(output_slice);
   buffer[1] = new Real_Slice(output_slice);
   Real_Slice *swap;  

   // --- COMPUTE PHANTOM SLICE WEIGHTING --- //

   double u, v, w0, w1, total_w = 0;

   if (z_start == z_stop){

      // --- SPECIAL CASE: sub-slice interpolation of phantom --- //
      // In this special case, both z_min and z_max lie within the same
      // phantom slice.  Only a single interpolation and interpolation
      // is required between phantom samples at the boundaries of the 
      // slice (i.e. at z_start and z_start+1).

      // Compute special case slice selection weighting
      u  = z_min - z_start; 
      v  = slice_thickness / phantom_thickness;
      w0 = v * (u + 0.5 * v);
      w1 = v - w0;
     
      // If slices do not exist in the phantom assume they are zero.
      if (z_start >= 0 && z_start < get_nslices()) {
         get_simulated_mag_phantom_slice(z_start, *(buffer[0]));
         output_slice.saxpy(w0, *(buffer[0]), output_slice); 
         total_w += w0;
      }
      if (z_start+1 >= 0 && z_start+1 < get_nslices()) { 
         get_simulated_mag_phantom_slice(z_start+1, *(buffer[1]));
         output_slice.saxpy(w1, *(buffer[1]), output_slice);
         total_w += w1;
      }
#ifdef TRACE
      cout << z_start << " " << w0 << " " << z_start+1 << " " << w1
           << " " << total_w << endl;
#endif   
   } else {

      // --- GENERAL CASE --- // 
      // In the general case at least two partial slices will exist.
      // One on each side of a phantom sample boundary.  Weights are
      // computed separately for each of these partial slices, and
      // then for any complete phantom slices between two phantom
      // samples.

      // --- FIRST PARTIAL SLICE --- //
      // The weighting of the first sample, w0, will be by the area
      // of the triangle with base u and height u (from interpolation).  
      // By similar triangles, the weight w1 will be the area of the
      // trapezoid with base u and heights u and 1.

      // Compute slice selection weighting for first partial slice
      u  = (z_start + 1) - z_min;
      w0 = 0.5 * u * u;
      w1 = u - w0;
   
      if (z_start >= 0 && z_start < get_nslices()) {
         get_simulated_mag_phantom_slice(z_start, *(buffer[0]));
         output_slice.saxpy(w0, *(buffer[0]), output_slice); 
         total_w += w0;
#ifdef TRACE
         cout << z_start << " " << w0 << endl;
#endif
      }
      if (z_start+1 >= 0 && z_start+1 < get_nslices()) { 
         get_simulated_mag_phantom_slice(z_start+1, *(buffer[1]));
         output_slice.saxpy(w1, *(buffer[1]), output_slice);
         total_w += w1;
#ifdef TRACE
         cout << z_start+1 << " " << w1 << endl;
#endif
      } else {
         buffer[1]->zeros();
      }

      // --- LOOP OVER COMPLETE SLICES --- //

      // Compute slice selection weighting for complete slices
      // i.e. as for first partial slice when u = 1.
      w0 = 0.5;  

      // Loop over all complete slices
      int z_n;
      for (z_n = z_start+1; z_n <= z_stop-1; z_n++){
   
         // Update buffers.
         // Buffer pointers are swapped so that the slice in buffer[1]
         // will be in buffer[0], and new slices are loaded in buffer[1].
         swap      = buffer[0];
         buffer[0] = buffer[1];
         buffer[1] = swap;
         if (z_n >= 0 && z_n < get_nslices()) {
            output_slice.saxpy(w0, *(buffer[0]), output_slice);
            total_w += w0;
#ifdef TRACE
            cout << z_n << " " << w0 << endl;
#endif
         }
         if (z_n+1 >= 0 && z_n+1 < get_nslices()) {
            get_simulated_mag_phantom_slice(z_n+1, *(buffer[1]));
            output_slice.saxpy(w0, *(buffer[1]), output_slice);
            total_w += w0;
#ifdef TRACE
            cout << z_n+1 << " " << w0 << endl;
#endif
         } else {
            buffer[1]->zeros();
         }
         
      }
   
      // --- LAST PARTIAL SLICE --- //
   
      // Compute slice selection weighting for last partial slice
      u  = z_max - z_stop;
      w1 = 0.5 * u * u;
      w0 = u - w1;

      // Update buffers
      swap      = buffer[0];
      buffer[0] = buffer[1];
      buffer[1] = swap;
      if (z_stop >= 0 && z_stop < get_nslices()) {
         output_slice.saxpy(w0, *(buffer[0]), output_slice);
         total_w += w0;
#ifdef TRACE
         cout << z_stop << " " << w0 << endl;
#endif
      }
      if (z_stop+1 >= 0 && z_stop+1 < get_nslices()) {
         get_simulated_mag_phantom_slice(z_stop+1, *(buffer[1]));
         output_slice.saxpy(w1, *(buffer[1]), output_slice); 
         total_w += w1;
#ifdef TRACE
         cout << z_stop+1 << " " << w1 << endl;
#endif
      }
  
   } 

   if (total_w > 0.0) {
      output_slice *= (1.0/total_w);
   }

#ifdef TRACE
   cout << "total_w: " << total_w << flush << endl;
#endif

   delete buffer[0];
   delete buffer[1];
}

//---------------------------------------------------------------------------
// Phantom::ideal_nn_slice_select
//
// Performs an ideal rectangular slice selection, for a given slice_thickness
// with the slice centred at z_centre.  Stores the result in output_slice.
// Phantom data is assumed to be nearest-neighbour interpolated between 
// samples and the slice selection profile is rectangular.
//
// Slice selection weighting is performed by integrating over the portion
// of the rectanngular interpolation kernel that falls within each phantom
// slice boundary.   
//---------------------------------------------------------------------------

void Phantom::ideal_nn_slice_select(double z_centre,
                                    double slice_thickness,
                                    Real_Slice& output_slice) {

#ifdef DEBUG
   // Check that output_slice is the right size
   assert(this->is_same_slice_size_as(output_slice));
   assert(slice_thickness > 0.0);
#endif

   // Get the slice separation thickness of the phantom.
   double phantom_thickness = get_voxel_step(SLICE);

#ifdef DEBUG
   assert(phantom_thickness > 0.0);
#endif

   // Compute slice coordinate limits
   // Phantom z-coordinate starts at 0.
   double z_min = (z_centre - 0.5*slice_thickness)/phantom_thickness;
   double z_max = (z_centre + 0.5*slice_thickness)/phantom_thickness;

   int z_start = (int)floor(z_min);
   int z_stop  = (int)floor(z_max);

   // z_start is rounded up at the half pixel
   // z_stop is rounded down at the half pixel
   if ((z_min - z_start) >= 0.5)
      z_start = z_start + 1;
   if ((z_max - z_stop) > 0.5)
      z_stop  = z_stop + 1;
   
   // Initialize phantom slice buffer
   output_slice.zeros();
   Real_Slice buffer(output_slice.get_nrows(), output_slice.get_ncols());

   // --- COMPUTE PHANTOM SLICE WEIGHTING --- //

   double w, total_w = 0;

   if (z_start == z_stop){

      // --- SPECIAL CASE: sub-slice interpolation of phantom --- //
      // In this special case, both z_min and z_max lie within a
      // phantom slice.  

      // Load simulated phantom slice into the buffer
      // If slices do not exist in the phantom assume they are zero.

      if (z_start >= 0 && z_start < get_nslices()) {

         get_simulated_mag_phantom_slice(z_start, buffer);

         // Compute special case slice selection weighting
         w = (z_max - z_min);
     
         // Weight phantom slice
         output_slice.saxpy(w, buffer, output_slice); 
         total_w += w;
      }
  
   } else {

      // --- GENERAL CASE ---//
      // In the general case at least two partial slices will exist.
      // One on each side of a phantom sample boundary.  Weights are
      // computed separately for each of these partial slices, and
      // then for any complete phantom slices between two phantom
      // samples.

      // --- FIRST PARTIAL SLICE --- //

      // Load simulated phantom slice into the buffer.
      // If slices do not exist in the phantom assume they are zero.
      if (z_start >= 0 && z_start < get_nslices()) {

         get_simulated_mag_phantom_slice(z_start, buffer);

         // Compute partial slice weight, and weight the slice.
         w = (z_start + 0.5 - z_min);
         output_slice.saxpy(w, buffer, output_slice);
         total_w += w;
      }
 
      // --- LOOP OVER COMPLETE SLICES --- //

      // Loop over all complete slices
      int z_n;
      for (z_n = z_start+1; z_n <= z_stop-1; z_n++){
     
         // Load simulated phantom slice into the buffer. 
         if (z_n >= 0 && z_n < get_nslices()) {

            get_simulated_mag_phantom_slice(z_n, buffer);

            // Weight the slice (w = 1)
            output_slice += buffer;
            total_w += 1;
         }

      }

      // --- LAST PARTIAL SLICE --- //

      // Load simulated phantom slice into the buffer.
      if (z_stop >= 0 && z_stop < get_nslices()) {

         get_simulated_mag_phantom_slice(z_stop, buffer);

         // Compute partial slice weight, and weight the slice.
         w = (z_max - z_stop + 0.5);
         output_slice.saxpy(w, buffer, output_slice);
         total_w += w;
      }

   }

   if (total_w > 0.0) { 
      output_slice *= (1.0/total_w);
   }

}

//---------------------------------------------------------------------------
// Phantom::compute_partial_volume
// Computes intra-slice partial volume by weighting phantom voxels
// assuming uniform intensity across the voxel.
//---------------------------------------------------------------------------

void Phantom::compute_partial_volume(const Real_Slice& phantom_slice,
                                     double row_step, double col_step,
                                     double row_shift, double col_shift,
                                     Real_Slice& image_slice) {


#ifdef DEBUG
   assert(row_step > 0.0);
   assert(col_step > 0.0);
#endif

   // --- Voxel geometry --- //

   const double phantom_row_step = get_voxel_step(ROW);
   const double phantom_col_step = get_voxel_step(COLUMN);

#ifdef DEBUG
   assert(phantom_row_step > 0.0);
   assert(phantom_col_step > 0.0);
#endif

   // --- Counters --- //

   double rmin, rmax, cmin, cmax;
   double rl, ru, cl, cu;
   int    irow, icol, prow, pcol;
   int    prow_start, prow_stop, pcol_start, pcol_stop;
   double data, row_weight;
   double total_weight;

   // --- LOOP OVER EACH VOXEL IN THE IMAGE --- //
   // Where possible, many computations are moved out of loops
   // for efficiency.

   // --- Initialize the row boundary variables --- //
   // Compute the initial values of rmin and rmax for the prow loop,
   // which give the boundaries of the current image row.

   rmin = row_shift - 0.5*row_step;
   rmax = rmin + row_step;

   // --- Loop over each row in the image --- //

   for (irow=0; irow<image_slice.get_nrows(); irow++){

      // --- Find all phantom rows within the current image row --- //
      // Row boundaries are used as limits for the prow loop

      prow_start = (int)floor((rmin+0.5*phantom_row_step)/phantom_row_step);
      prow_stop  = (int) ceil((rmax+0.5*phantom_row_step)/phantom_row_step);

      // --- Initialize the column boundary variables --- //
      // Compute the initial values of cmin and cmax for the pcol loop,
      // which give the boundaries of the current image column.

      cmin  = col_shift - 0.5*col_step;
      cmax  = cmin + col_step;

      // --- Loop over each column in the image --- //

      for (icol=0; icol<image_slice.get_ncols(); icol++){

         // --- FOR EACH VOXEL IN THE IMAGE --- //

         data = 0;
         total_weight = 0;

#ifdef DEBUG
         // Loop invariant (relaxed for fp ops)
         assert(rmax - rmin - row_step < 1E-5);
         assert(cmax - cmin - col_step < 1E-5);
#endif

         // --- Find all phantom columns within the current image column --- //
         // Column boundaries are used as the limits for the pcol loop.

         pcol_start = (int)floor((cmin+0.5*phantom_col_step)/phantom_col_step);
         pcol_stop  = (int) ceil((cmax+0.5*phantom_col_step)/phantom_col_step);

         // --- COMPUTE VOXEL WEIGHTING --- //

         // Loop over each phantom voxel which overlaps with
         // the current image voxel.

         // --- Initialize the phantom row boundary variables --- //
         // ru and rl give the upper and lower boundaries of a phantom voxel
         // used in the prow loop.

         rl    = rmin;
         ru    = (prow_start + 0.5)*phantom_row_step;
         if (ru > rmax) ru = rmax;

         for (prow=prow_start; prow<prow_stop; prow++) {
             
#ifdef DEBUG
            // Loop invariant (relaxed for fp ops)
            assert(ru >= rl);
            assert(ru - rl - phantom_row_step <= 1E-5);
#endif
 
            row_weight = ru - rl;

            // --- Initialize the phantom column boundary variables --- //
            // cl and cu give the lower and upper boundaries of a phantom voxel
            // used in the pcol loop.

            cl    = cmin;
            cu    = (pcol_start + 0.5)*phantom_col_step;
            if (cu > cmax) cu = cmax;

            for (pcol=pcol_start; pcol<pcol_stop; pcol++) {

#ifdef DEBUG
               // Loop invariant (relaxed for fp ops)
               assert(cu >= cl);
               assert(cu - cl - phantom_col_step <= 1E-5);
#endif

               if (prow >= 0 && prow < get_nrows() &&
                   pcol >= 0 && pcol < get_ncols() ) {
                  data += row_weight * (cu - cl) * phantom_slice(prow, pcol);
                  total_weight += row_weight * (cu - cl);
               }

               // Update phantom column boundaries
 
               cl = cu;
               cu += phantom_col_step;
               if (cu > cmax) cu = cmax;

            }

            // Update phantom row boundaries

            rl = ru;
            ru += phantom_row_step; 
            if (ru > rmax) ru = rmax;

         } 

         // --- SAVE THE WEIGHTED VOXEL --- //

         //image_slice.set_value(irow, icol, data/total_weight);
         if (total_weight > 0.0) {
            image_slice(irow, icol) = data/total_weight;
         } else {
            image_slice(irow, icol) = 0.0;
         }

#ifdef DEBUG
         // Check weight computation
         assert(total_weight - (row_step*col_step) <= 1E-5);
         assert(total_weight >= 0.0);
#endif

         // --- Update the column boundary variables --- //

         cmin  = cmax;
         cmax += col_step; 

      } 

      // --- Update the row boundary variables --- //

      rmin  = rmax;
      rmax += row_step;

   } 

}

//---------------------------------------------------------------------------
// Phantom::_direct_fourier_resample
// Resample the simulated slice's Fourier transform using direct evaluation.
//---------------------------------------------------------------------------

void Phantom::_direct_fourier_resample(const Complex_Slice& /*sim_slice*/,
                                       Complex_Slice& /*raw_slice*/) {

#ifdef DEBUG
   // assert that raw_slice is the right size
#endif

   cerr << "Phantom::_direct_fourier_resample not implemented." << endl;

}

//---------------------------------------------------------------------------
// Phantom::_fft_fourier_resample
// Resample the simulated slice's Fourier transform using FFT evaluation.
// Works only for power of two matrix sizes.
// Does not work when image bandwidth is larger than the sampled phantom's
// bandwidth.  In that case use Phantom::_direct_fourier_resample.
//---------------------------------------------------------------------------

void Phantom::_fft_fourier_resample(const Complex_Slice& /*sim_slice*/,
                                    Complex_Slice& /*raw_slice*/) {

#ifdef DEBUG
   // assert that raw_slice is the right size
#endif

   cerr << "Phantom::_fft_fourier_resample not implemented." << endl;

}

//---------------------------------------------------------------------------
// Phantom::initialize_chirp
// Precomputes Chirp DFT filters required for the row-column decomposition
// of the 2-D Chirp DFT.
//---------------------------------------------------------------------------

void Phantom::initialize_chirp(unsigned int out_row_length,
                               unsigned int out_col_length,
                               double out_row_fov,
                               double out_col_fov) {

   // The 2-D Chirp DFT is performed by a row-column decomposition of
   // 1-D Chirp DFTs.  A 1-D Chirp is performed on each row of the matrix
   // followed by a 1-D Chirp on each column of the matrix.
   // Note:   a row is a 1-D slice taken along the COLUMN direction, and
   // a col is a 1-D slice taken along the ROW direction.   There are
   // thus get_nrows() rows each of length get_ncols() elements and 
   // intersample distance get_voxel_step(COLUMN).

   const unsigned int in_row_length = get_ncols();
   const unsigned int in_col_length = get_nrows();
   const double in_row_fov    = in_row_length * get_voxel_step(COLUMN);
   const double in_col_fov    = in_col_length * get_voxel_step(ROW);

   // Compute chirp DFT frequency parameters
   // The input (phantom) sampling determines the Nyquist frequency and
   // scale of the frequency space, all frequency parameters must be given
   // in terms of sampled frequency units (in the range -M_PI..M_PI).  
   // The continuous freq separation of the input is k_i = 1/FOV_i, while   
   // the sampled freq separation is k_is = 2*M_PI/N_i.   
   // The conversion from continuous to sampled frequencies to is then
   // 2*M_PI*FOV_i/N_i;

   const double row_w_step    = 2*M_PI * in_row_fov / out_row_fov /
                                in_row_length;

   const double col_w_step    = 2*M_PI * in_col_fov / out_col_fov /
                                in_col_length;

   // The frequency samples are placed at -(N_o/2)..(N_o/2-1).  Multiplying
   // by the frequency separation (w_step) gives the position of the initial
   // frequency sample.

   const double row_w_initial = -row_w_step*(double)out_row_length/2.0;
   const double col_w_initial = -col_w_step*(double)out_row_length/2.0;

   // Precompute the Chirp DFT filters needed to compute 1-D transforms
   // over row or col slices.   The 2-D Chirp DFT is formed by a
   // row-column decomposition using the 1-D DFTs.

   if (row_chirp != NULL) delete row_chirp;
   row_chirp = new Chirp_Algorithm(in_row_length, out_row_length,
                                   row_w_initial, row_w_step);

   if (col_chirp != NULL) delete col_chirp;
   col_chirp = new Chirp_Algorithm(in_col_length, out_col_length,
                                   col_w_initial, col_w_step);

   if (row_weight != NULL) delete row_weight;
   row_weight = new double[2*out_row_length];
   _compute_chirp_weight(row_weight, in_row_length, out_row_length, 
                         in_row_fov, out_row_fov);

   if (col_weight != NULL) delete col_weight;
   col_weight = new double[2*out_col_length];
   _compute_chirp_weight(col_weight, in_col_length, out_col_length, 
                         in_col_fov, out_col_fov);

   // A temporary slice to store the partial result of the row-wise 
   // 1-D chirp DFT

   if (tmp_slice != NULL) delete tmp_slice;
   tmp_slice = new Complex_Slice(in_row_length, out_col_length);

}

//---------------------------------------------------------------------------
// Phantom::_chirp_fourier_resample 
// Resample the simulated slice's Fourier transform using a Chirp DFT.
// Works for any arbitrary matrix size.
// Does not work when image bandwidth is larger than the sampled phantom's
// bandwidth.  In that case use Phantom::_direct_fourier_resample.
//---------------------------------------------------------------------------

void Phantom::_chirp_fourier_resample(const Real_Slice& sim_slice,
                                      Complex_Slice& raw_slice) {

#ifdef DEBUG
   // Ensure that chirps have been initialized
   assert(row_chirp != NULL);
   assert(col_chirp != NULL);
   assert(tmp_slice != NULL);

   // assert that sim_slice is the right size
   assert(sim_slice.get_ncols() == row_chirp->get_input_length());
   assert(sim_slice.get_nrows() == col_chirp->get_input_length());

   // assert that raw_slice is the right size
   assert(raw_slice.get_ncols() == row_chirp->get_output_length());
   assert(raw_slice.get_nrows() == col_chirp->get_output_length());
#endif

   unsigned int n;
   for (n=0; n<sim_slice.get_nrows(); n++){
      row_chirp->apply(FALSE, sim_slice.row_ptr(n), sim_slice.get_row_stride(),
                       tmp_slice->row_ptr(n), tmp_slice->get_row_stride());
//      _apply_weight(tmp_slice->row_ptr(n), tmp_slice->get_row_stride(),
//                    row_weight, tmp_slice->get_row_length());
   }
   for (n=0; n<raw_slice.get_ncols(); n++){
      col_chirp->apply(TRUE, tmp_slice->col_ptr(n), tmp_slice->get_col_stride(),
                       raw_slice.col_ptr(n), raw_slice.get_col_stride());
//      _apply_weight(raw_slice.col_ptr(n), raw_slice.get_col_stride(),
//                    col_weight, raw_slice.get_col_length());
   }

   for (n=0; n<raw_slice.get_nrows(); n++){
      _apply_weight(raw_slice.row_ptr(n), raw_slice.get_row_stride(),
                    row_weight, raw_slice.get_row_length());
   }
   for (n=0; n<raw_slice.get_ncols(); n++){
      _apply_weight(raw_slice.col_ptr(n), raw_slice.get_col_stride(),
                    col_weight, raw_slice.get_col_length());
   }

}

//---------------------------------------------------------------------------
// Phantom::_compute_chirp_weight
//---------------------------------------------------------------------------

void Phantom::_compute_chirp_weight(double weight[], 
                                    unsigned int in_length, 
                                    unsigned int out_length,
                                    double in_fov, double out_fov) {

   const double shift      = M_PI_2*(1.0 - (in_fov/out_fov)); 
   const double voxel_size = in_fov/(double)in_length;
   const double factor     = 0.5*voxel_size/out_fov;
 
   double angle, linear_comp;

   int n;
   for (n=0; n<2*out_length; n+=2){
      angle = shift*n;
      linear_comp =  voxel_size * 
                     _sinc(factor*(n-(double)out_length)) *
                     _sinc(factor*(n-(double)out_length));
      weight[n]   =  linear_comp*cos(angle);
      weight[n+1] = -linear_comp*sin(angle);
   }

}

//---------------------------------------------------------------------------
// Phantom::_apply_weight
//---------------------------------------------------------------------------

void Phantom::_apply_weight(Real_Scalar in[], unsigned int stride,
                            double weight[], unsigned int nelements) {

   int m, n;
   double tempr, tempi, filtr, filti;

   for(n=0, m=0; n<2*nelements; n+=2, m+=(2*stride)){
      tempr = in[m];
      tempi = in[m+1];
      filtr = weight[n];
      filti = weight[n+1];

      in[m]   = (Real_Scalar)(tempr * filtr - tempi * filti);
      in[m+1] = (Real_Scalar)(tempr * filti + tempi * filtr);
   }
}

//---------------------------------------------------------------------------
// Phantom::_compensate_for_linear_kernel
// Weights Fourier samples to compensate for the assumed underlying
// linear interpolation kernel.
//---------------------------------------------------------------------------

void Phantom::_compensate_for_linear_kernel(Complex_Slice& raw_slice) {
 
   unsigned int m, n;
   double       kr, kc, weightr, weightc;

   const double row_step = get_voxel_step(ROW);
   const double col_step = get_voxel_step(COLUMN);
   const double kr_step  = 1.0/_pseq->get_fov(ROW);
   const double kc_step  = 1.0/_pseq->get_fov(COLUMN);
   const double kr0      = -((double)raw_slice.get_nrows()*kr_step/2.0);
   const double kc0      = -((double)raw_slice.get_ncols()*kc_step/2.0);
 
   for (m=0, kr=kr0; m<raw_slice.get_nrows(); m++, kr+=kr_step){
      weightr = row_step * _sinc(row_step*kr) * _sinc(row_step*kr);
      for (n=0, kc=kc0; n<raw_slice.get_ncols(); n++, kc+=kc_step){
         weightc = weightr * col_step * 
                  _sinc(col_step*kc) * _sinc(col_step*kc);

         raw_slice.real(m,n) *= weightc;
         raw_slice.imag(m,n) *= weightc;
      }
   }

}
