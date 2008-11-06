//==========================================================================
// IMINCFILE.CXX
// Member functions for class I_MINC_File.
// Inherits from:  MINC_File
// Base class to:  IO_MINC_File
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/imincfile.cxx,v 1.2 2008-11-06 10:58:22 rotor Exp $
// $Log: imincfile.cxx,v $
// Revision 1.2  2008-11-06 10:58:22  rotor
//  * fixed includes for iostream and friends
//  * updated for new release (1.0.2)
//
// Revision 1.1  2003/05/30 16:43:08  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:26:48  rkwan
// Release 2.5
//
// Revision 1.2  1995/12/11  14:15:31  rkwan
// Updated for fast_iso_model.
//
// Revision 1.1  1995/10/13  14:24:55  rkwan
// Initial revision
//
//==========================================================================

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream>

#include "imincfile.h"
#include "mriimage.h"

//--------------------------------------------------------------------------
// I_MINC_File constructors
//--------------------------------------------------------------------------

I_MINC_File::I_MINC_File() : MINC_File() {}

//--------------------------------------------------------------------------
// I_MINC_File destructor 
//--------------------------------------------------------------------------

I_MINC_File::~I_MINC_File() {}

//--------------------------------------------------------------------------
// I_MINC_File::open
// Opens a MINC volume for reading.  Sets status bits on failure, 
// and returns the id of the open MINC file.
//--------------------------------------------------------------------------

int I_MINC_File::open(const char *path) {

   // If file is already open, close it
   if (this->is_open()){
      this->close();
   }

   // Open the file
   set_error_off();
   if ((_MINCid = miopen((char *)path, NC_NOWRITE)) == MI_ERROR){
      _set_failure();              // Could not open file
   } else {
      _setstate(file_open_bit);    // Otherwise, file is open
      strcpy(_filepath, path);     // so save file name
   }
   set_error_on();

   return _MINCid;

}
  
//--------------------------------------------------------------------------
// I_MINC_File::attach_icv
// Attachs the icv to the MINC volume after all icv properties have
// been modified.
//--------------------------------------------------------------------------

int I_MINC_File::attach_icv(void) {
   int status;

   set_error_off();
   if ((status = miicv_attach(_icvid, _MINCid, ncvarid(_MINCid, MIimage))) 
       == MI_ERROR){
      _set_failure();                 // could not attach icv
   } else {
      _setstate(icv_attached_bit);    // icv is attached
      _get_volume_info();             // load volume information from disk
   }
   set_error_on();

   return status;
}

//--------------------------------------------------------------------------
// I_MINC_File::load_hyperslab
// Reads a hyperslab of data into memory.   start specifies the beginning
// coordinates of the hyperslab, and count specifies the number of
// elements in each dimension of the hyperslab.
//--------------------------------------------------------------------------

int I_MINC_File::load_hyperslab(long start[], long count[], void *volume){

#ifdef DEBUG
   assert(this->is_good());
   assert(volume != NULL);
#endif
   return miicv_get(_icvid, start, count, volume);

}

//--------------------------------------------------------------------------
// I_MINC_File::load_slice
// Reads a slice of data into memory. 
//--------------------------------------------------------------------------

int I_MINC_File::load_slice(int slice_num, void *slice){

   long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS];
   int  ndims;

   // Get number of dimensions
   ndims = _volume_info.number_of_dimensions;

   // Set up the start and count variables for reading the volume
   (void) miset_coords(MAX_VAR_DIMS, 0, start);

   if (ndims >= 3){
      start[ndims-3] = slice_num;
      count[ndims-3] = 1;
   }
   count[ndims-2] = _volume_info.length[ndims-2];
   count[ndims-1] = _volume_info.length[ndims-1];

   // Read in the volume
   return this->load_hyperslab(start,count,slice);

}

int I_MINC_File::load_slice(int slice_num, MRI_Image& image){

   long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS];
   int  ndims;
   double min, max;

#ifdef DEBUG
   assert(this->is_good());
#endif

   // Get number of dimensions
   ndims = _volume_info.number_of_dimensions;

   // Set up the start and count variables for reading the volume
   (void) miset_coords(3, 0, start);

   if (ndims >= 3){
      start[ndims-3] = slice_num;
      count[ndims-3] = 1;
   }
   count[ndims-2] = _volume_info.length[ndims-2];
   count[ndims-1] = _volume_info.length[ndims-1];

   // Read in the image maximum and minimum
   if (mivarget1(_MINCid, ncvarid(_MINCid, MIimagemin), start, NC_DOUBLE, 
                 NULL, &min) == MI_ERROR){
      min = 0;
      cerr << "WARNING: Failed to read image-min in I_MINC_File::load_slice."
           << endl 
           << "         Assuming image-min == 0." << flush << endl;
   } 
   if (mivarget1(_MINCid, ncvarid(_MINCid, MIimagemax), start, NC_DOUBLE, 
                 NULL, &max) == MI_ERROR){
      max = 4095;
      cerr << "WARNING: Failed to read image-max in I_MINC_File::load_slice."
           << endl 
           << "         Assuming image-max == 4095." << flush << endl;
   }

   image.set_real_minimum(min);
   image.set_real_maximum(max);

   // Read in the volume
   return this->load_hyperslab(start,count,(void *)image);

}

//--------------------------------------------------------------------------
// I_MINC_File protected member functions
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// I_MINC_File::_get_volume_info
// Reads volume information from the MINC file into the internal
// information structure. 
//--------------------------------------------------------------------------

void I_MINC_File::_get_volume_info(void){
   int imgid, varid;
   int idim, ndims;
   int dim[MAX_VAR_DIMS];
   char *dimname;
   int att_length, sign_index;

   if (!this->is_good()) {
      cerr << "I_MINC_File::_get_volume_info: "
           << "File must be open and have an attached ICV." << flush << endl;
      exit(EXIT_FAILURE);
   }

   // Get the minc file id and the image variable id

   (void) miicv_inqint(_icvid, MI_ICV_VARID, &imgid);
   if ((_MINCid == MI_ERROR) || (imgid == MI_ERROR)) {
      cerr << "Could not read image variable." << flush << endl;
      exit(EXIT_FAILURE);
   }

   // Get the volume min and max 
//   (void) miicv_inqdbl(_icvid, MI_ICV_NORM_MIN, &_volume_info.minimum);
//   (void) miicv_inqdbl(_icvid, MI_ICV_NORM_MAX, &_volume_info.maximum);

   // Get the list of dimensions subscripting the image variable 
   (void) ncvarinq(_MINCid, imgid, NULL, NULL, &ndims, dim, NULL);
   (void) miicv_inqint(_icvid, MI_ICV_NUM_DIMS, &ndims);
   _volume_info.number_of_dimensions = ndims;

   // Initialize step and start to unknown values
   for (idim=0; idim<ndims; idim++){
      _volume_info.step[idim] =  HUGE_VAL;
      _volume_info.start[idim] = HUGE_VAL;
   }

   // Turn off error processing
   set_error_off();

   // Loop through the dimensions checking them and getting their sizes
   for (idim=0; idim<ndims; idim++){
      
      // Get pointers to the appropriate dimension size and name
      dimname = _volume_info.dimension_names[idim];
      (void)ncdiminq(_MINCid,dim[idim],dimname,&(_volume_info.length[idim]));

      // Get step and start directly from the variables
      if ((varid = ncvarid(_MINCid, _volume_info.dimension_names[idim]))
          != MI_ERROR){
         (void)miattget1(_MINCid, varid, MIstep, NC_DOUBLE, 
                         &(_volume_info.step[idim]));
         (void)miattget1(_MINCid, varid, MIstart, NC_DOUBLE,
                         &(_volume_info.start[idim]));
      } 
   }

   // Get data type
   (void)ncvarinq(_MINCid, imgid, NULL, &_volume_info.datatype,
                  &ndims, dim, NULL);

   // Get signtype
   if ((miattgetstr(_MINCid, imgid, MIsigntype, MI_MAX_ATTSTR_LEN, 
                    _volume_info.signtype) == NULL) ||
       ((strcmp(_volume_info.signtype, MI_UNSIGNED)!=0) &&
        (strcmp(_volume_info.signtype, MI_SIGNED)!=0))){
      if (_volume_info.datatype == NC_BYTE)
         (void) strcpy(_volume_info.signtype, MI_UNSIGNED);
      else
         (void) strcpy(_volume_info.signtype, MI_SIGNED);
   }
   sign_index = (strcmp(_volume_info.signtype, MI_UNSIGNED)==0) ? 0 : 1;

   // Get valid range
   if ((miattget(_MINCid, imgid, MIvalid_range, NC_DOUBLE, 2,
                 _volume_info.valid_range, &att_length) == MI_ERROR) ||
       (att_length != 2)){
      if (miattget1(_MINCid, imgid, MIvalid_min, NC_DOUBLE,
          &_volume_info.valid_range[0]) == MI_ERROR)
         _volume_info.valid_range[0] = 
            MINC_File::default_min[_volume_info.datatype][sign_index];
      if (miattget1(_MINCid, imgid, MIvalid_max, NC_DOUBLE,
          &_volume_info.valid_range[1]) == MI_ERROR)
         _volume_info.valid_range[1] =
            MINC_File::default_max[_volume_info.datatype][sign_index];
   }

   set_error_on();

}






