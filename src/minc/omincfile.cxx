//==========================================================================
// OMINCFILE.CXX
// Member functions for class O_MINC_File.
// Inherits from:  MINC_File
// Base class to:  IO_MINC_File
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/omincfile.cxx,v 1.1 2003-05-30 16:43:09 bert Exp $
// $Log: omincfile.cxx,v $
// Revision 1.1  2003-05-30 16:43:09  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.6  1996/05/29  19:05:10  rkwan
// GNU warning fix.
//
// Revision 2.5  1996/05/29  16:29:30  rkwan
// Release 2.5
//
// Revision 1.3  1996/01/17  17:59:49  rkwan
// Update for tx/rx coil modelling.
//
// Revision 1.2  1995/12/11  14:18:36  rkwan
// Updated for fast_iso_model.
//
//==========================================================================

#include "imincfile.h"
#include "omincfile.h"
#include <string.h>

//--------------------------------------------------------------------------
// O_MINC_File constructor 
//--------------------------------------------------------------------------

O_MINC_File::O_MINC_File() : MINC_File() {}

//--------------------------------------------------------------------------
// O_MINC_File destructor
//--------------------------------------------------------------------------

O_MINC_File::~O_MINC_File() {}

//--------------------------------------------------------------------------
// O_MINC_File::open
// Opens a MINC volume for writing.  Sets status bits on failure, 
// and returns the id of the open MINC file.
//--------------------------------------------------------------------------

int O_MINC_File::open(const char *path) {

   // If file is already open close it
   if (this->is_open()){
      this->close();
   }

   // Open the file
   set_error_off();
   if ((_MINCid = miopen((char *)path, NC_WRITE)) == MI_ERROR){
      _set_failure();                // Could not open file
   } else {
      _setstate(file_open_bit);      // Otherwise, file is open
      strcpy(_filepath, path);       // so save file name
   }
   set_error_on();

   return _MINCid;

}

//--------------------------------------------------------------------------
// O_MINC_File::create
// Creates a new MINC file for writing.
//--------------------------------------------------------------------------

int O_MINC_File::create(const char *path, int mode) {

   // If file is already open close it
   if (this->is_open()){
      this->close();
   }

   // Open the file
   set_error_on();
   if ((_MINCid = micreate((char *)path, mode)) == MI_ERROR){
      _set_failure();               // Could not create file
   } else {
      _setstate(file_open_bit);     // Otherwise, file has been created
      strcpy(_filepath, path);      // so save path name
   }
   set_error_off();

   return _MINCid;

}

//--------------------------------------------------------------------------
// O_MINC_File::save_hyperslab
// Writes a hyperslab of data to the MINC file.  start specifies the
// beginning coordinates of the hyperslab, and count specifies the 
// number of elements in each dimension of the hyperslab.
//--------------------------------------------------------------------------

int O_MINC_File::save_hyperslab(long start[], long count[], void *volume){

#ifdef DEBUG
   assert(this->is_good());
   assert(volume != NULL);
#endif

   return miicv_put(_icvid, start, count, volume);

}

//--------------------------------------------------------------------------
// O_MINC_File::save_slice
// Writes a slice of data to the MINC file.
//--------------------------------------------------------------------------

int O_MINC_File::save_slice(int slice_num, void *slice){

   long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS];
   int  ndims = _volume_info.number_of_dimensions;

#ifdef DEBUG
   assert(this->is_good());
   assert(slice != NULL);
#endif

   // Set up the start and count variables
   (void) miset_coords(MAX_VAR_DIMS, 0, start);

   if (ndims >= 3){
      start[ndims-3]  = slice_num;
      count[ndims-3]  = 1;
   }
   count[ndims-2] = _volume_info.length[ndims-2];
   count[ndims-1] = _volume_info.length[ndims-1];
   
   // Write out slice min and max
   (void) mivarput1(_MINCid, ncvarid(_MINCid, MIimagemin), start,
                    NC_DOUBLE, NULL, &_volume_info.valid_range[0]);
   (void) mivarput1(_MINCid, ncvarid(_MINCid, MIimagemax), start,
                    NC_DOUBLE, NULL, &_volume_info.valid_range[1]);

   // Write out the slice 
   return this->save_hyperslab(start, count, slice);

}

int O_MINC_File::save_slice(int slice_num, MRI_Image& image){

   long start[MAX_VAR_DIMS], count[MAX_VAR_DIMS];
   int  ndims = _volume_info.number_of_dimensions;
   double slice_min, slice_max;

#ifdef DEBUG
   assert(this->is_good());
#endif

   // Set up the start and count variables
   (void) miset_coords(MAX_VAR_DIMS, 0, start);
   if (ndims >= 3){ 
      start[ndims-3]  = slice_num;
      count[ndims-3]  = 1;
   }
   count[ndims-2] = image.get_nrows();
   count[ndims-1] = image.get_ncols();

   // Write out slice min and max
   slice_min = image.get_real_minimum();
   slice_max = image.get_real_maximum();

   (void) mivarput1(_MINCid, ncvarid(_MINCid, MIimagemin), start,
                    NC_DOUBLE, NULL, &slice_min);
   (void) mivarput1(_MINCid, ncvarid(_MINCid, MIimagemax), start,
                    NC_DOUBLE, NULL, &slice_max);

   // Write out the slice
   return this->save_hyperslab(start, count, (void *)image);

}

//--------------------------------------------------------------------------
// O_MINC_File::set_volume_info
// Updates the MINC file's volume information.
//--------------------------------------------------------------------------

void O_MINC_File::set_volume_info(const I_MINC_File& ifile, 
                                  const char *arg_string){

   if (ifile.is_bad()) {
      cerr << "O_MINC_File::set_volume_info: bad input MINC file specified."
           << flush << endl;
      exit(EXIT_FAILURE);
   }
   this->set_volume_info(ifile.get_cdfid(), ifile._volume_info, arg_string);

}

void O_MINC_File::set_volume_info(int inMINCid,
                                  const Volume_Info &volume_info,
                                  const char *arg_string){

   int dim[MAX_VAR_DIMS], idim, varid;
   int excluded_vars[10], nexcluded;
   int ndims = volume_info.number_of_dimensions;

   // Check if the MINC volume is bad
   if (this->is_bad()){
      cerr << "O_MINC_File::set_volume_info: bad output MINC file specified."
           << flush << endl;
      exit(EXIT_FAILURE);
   }
   
   // --- Copy volume info --- //

//   _volume_info.minimum = volume_info.minimum;
//   _volume_info.maximum = volume_info.maximum;

   _volume_info.number_of_dimensions = volume_info.number_of_dimensions;
   _volume_info.datatype             = volume_info.datatype;
   _volume_info.valid_range[0]       = volume_info.valid_range[0];
   _volume_info.valid_range[1]       = volume_info.valid_range[1];

   strcpy(_volume_info.signtype, volume_info.signtype);

   for(idim = 0; idim < ndims; idim++){
      _volume_info.length[idim] = volume_info.length[idim];
      _volume_info.step[idim]   = volume_info.step[idim];
      _volume_info.start[idim]  = volume_info.start[idim];

      strcpy(_volume_info.dimension_names[idim], 
              volume_info.dimension_names[idim]);
   }

   // --- Create the dimensions --- //

   for(idim = 0; idim < ndims; idim++){
      dim[idim] = ncdimdef(_MINCid, volume_info.dimension_names[idim],
                           volume_info.length[idim]);
   }

   // If an input file is provided, copy all header info from that file except
   // image, image-max, image-min
   if (inMINCid != MI_ERROR) {

      // Look for the image variable and the image-max/min variables so that
      // we can exclude them from the copy
      nexcluded = 0;
      excluded_vars[nexcluded] = ncvarid(inMINCid, MIimage);
      if (excluded_vars[nexcluded] != MI_ERROR) nexcluded++;
      excluded_vars[nexcluded] = ncvarid(inMINCid, MIimagemax);
      if (excluded_vars[nexcluded] != MI_ERROR) nexcluded++;
      excluded_vars[nexcluded] = ncvarid(inMINCid, MIimagemin);
      if (excluded_vars[nexcluded] != MI_ERROR) nexcluded++;

      // Copy the variable definitions
      (void) micopy_all_var_defs(inMINCid, _MINCid, nexcluded, excluded_vars);

   }

   // Set up the dimension variables. If the variable doesn't exist, create
   // it (either no input file or variable did not exist in it). If the
   // dimensions are not standard, then no variable is created.

   for (idim=0; idim < ndims; idim++) {
      set_error_off();
      varid = ncvarid(_MINCid, volume_info.dimension_names[idim]);
      if (varid == MI_ERROR) {
         varid = micreate_std_variable(_MINCid,
                                    (char *)volume_info.dimension_names[idim],
                                    NC_LONG, 0, NULL);
      }
      set_error_on();

      if (varid != MI_ERROR) {
         (void) miattputdbl(_MINCid, varid, MIstep,
                            volume_info.step[idim]);
         (void) miattputdbl(_MINCid, varid, MIstart,
                            volume_info.start[idim]);
      }
   }

   // Create the image, image-max and image-min variables
   _setup_image_variables(inMINCid, ndims, dim);

   // Add the time stamp to the history
   if (arg_string != NULL){
      _update_history(arg_string);
   }

   // Put the file in data mode
   (void) ncendef(_MINCid);

   // Copy over variable values
   if (inMINCid != MI_ERROR) {
      (void) micopy_all_var_values(inMINCid, _MINCid,
                                   nexcluded, excluded_vars);
   }

   
}

//--------------------------------------------------------------------------
// O_MINC_File Protected member functions
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// O_MINC_File::_setup_image_variables
// Updates image variables in the MINC file.
//--------------------------------------------------------------------------

void O_MINC_File::_setup_image_variables(int inMINCid, int ndims, int dim[]){

   int imgid, maxid, minid;

   // Create the image variable, copy attributes, set the signtype attribute,
   // set the valid range attribute and delete valid max/min attributes
   imgid = micreate_std_variable(_MINCid, MIimage, _volume_info.datatype, 
                                 ndims, dim);
   if (inMINCid != MI_ERROR) {
      (void) micopy_all_atts(inMINCid, ncvarid(inMINCid, MIimage),
                        _MINCid, imgid);
      (void) ncattdel(_MINCid, imgid, MIvalid_max);
      (void) ncattdel(_MINCid, imgid, MIvalid_min);
   }
   (void) miattputstr(_MINCid, imgid, MIsigntype, _volume_info.signtype);
   (void) ncattput(_MINCid, imgid, MIvalid_range, 
                   NC_DOUBLE, 2, _volume_info.valid_range);

   // Create the image max and min variables (varying over slices)
   maxid = micreate_std_variable(_MINCid, MIimagemax, NC_DOUBLE, 1, dim);
   minid = micreate_std_variable(_MINCid, MIimagemin, NC_DOUBLE, 1, dim);
   if (inMINCid != MI_ERROR) {
      (void) micopy_all_atts(inMINCid, ncvarid(inMINCid, MIimagemax),
                             _MINCid, maxid);
      (void) micopy_all_atts(inMINCid, ncvarid(inMINCid, MIimagemin),
                             _MINCid, minid);
   }

}

//--------------------------------------------------------------------------
// O_MINC_File::_update_history
// Update the history time stamp of the MINC file.
//--------------------------------------------------------------------------

void O_MINC_File::_update_history(const char *arg_string) {

   nc_type datatype;
   int att_length;
   char *string;

   // Get the history attribute length
   set_error_off();
   if ((ncattinq(_MINCid, NC_GLOBAL, MIhistory, &datatype, &att_length) 
        == MI_ERROR) || (datatype != NC_CHAR)) {
      att_length = 0;
   }
   set_error_on();
   att_length += strlen(arg_string) + 1;

   // Allocate a string and get the old history
   string = new char[att_length];
   string[0] = '\0';
   (void) miattgetstr(_MINCid, NC_GLOBAL, MIhistory, att_length,
                      string);

   // Add the new command and put the new history.
   (void) strcat(string, arg_string);
   (void) miattputstr(_MINCid, NC_GLOBAL, MIhistory, string);
   delete[] string;

}
