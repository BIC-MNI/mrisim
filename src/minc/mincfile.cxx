//==========================================================================
// MINCFILE.CXX
// Base class for MINC file classes.
// Inherits from:  MINC_ICV
// Base class to:  I_MINC_File, O_MINC_File
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/mincfile.cxx,v 1.2 2004-08-10 15:20:17 bert Exp $
// $Log: mincfile.cxx,v $
// Revision 1.2  2004-08-10 15:20:17  bert
// Comment out obsolete NC_xxx types
//
// Revision 1.1  2003/05/30 16:43:09  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:22:15  rkwan
// Release 2.5
//
// Revision 1.2  1995/12/11  14:16:14  rkwan
// Updated for fast_iso_model.
//
//==========================================================================

#include <string.h>
#include <limits.h>
#include <math.h>
#include "mincfile.h"

//--------------------------------------------------------------------------
// MINC_File::_ncoldopts
// Stores error checking option.
//--------------------------------------------------------------------------

int MINC_File::_ncoldopts =get_ncopts();

//--------------------------------------------------------------------------
// MINC_File::type 
// static look up table for data type names.
//--------------------------------------------------------------------------

char *MINC_File::type[] = {
   NULL, "byte", "char", "short", "long", "float", "double" 
};

//--------------------------------------------------------------------------
// MINC_File::default_min
// static look up table for minimum valid range.
//--------------------------------------------------------------------------

double MINC_File::default_min[][2] = {
      0.0, 0.0,
      0.0, (double)SCHAR_MIN,
      0.0, 0.0,
      0.0, (double)SHRT_MIN,
      0.0, (double)LONG_MIN,
      0.0, 0.0,
      0.0, 0.0 };

//--------------------------------------------------------------------------
// MINC_File::default_max
// static look up table for maximum valid range.
//--------------------------------------------------------------------------

double MINC_File::default_max[][2] = {
      0.0, 0.0,
      (double)UCHAR_MAX, (double)SCHAR_MAX,
      0.0, 0.0,
      (double)USHRT_MAX, (double)SHRT_MAX,
      (double)ULONG_MAX, (double)LONG_MAX,
      1.0, 1.0,
      1.0, 1.0 };

//--------------------------------------------------------------------------
// MINC_File constructors
//--------------------------------------------------------------------------

MINC_File::MINC_File() : MINC_ICV() {
   _state = 0;
}

//--------------------------------------------------------------------------
// MINC_File destructor
//--------------------------------------------------------------------------

MINC_File::~MINC_File() {

   if (this->is_open()){
      this->close();
   }
   if (this->has_ICV_attached()){
      this->detach_icv();
   }

}

//--------------------------------------------------------------------------
// MINC_File::open
// Opens a MINC volume.   Sets status bits on failure, and returns
// the id of the open MINC file.
//--------------------------------------------------------------------------

int MINC_File::open(const char *path, int mode) {

   // If file is already open close it
   if (this->is_open()){
      this->close();
   }

   // Open the file
   set_error_off();
   if ((_MINCid = miopen((char *)path, mode)) == MI_ERROR){
      _set_failure();             // Signal could not open file
   } else {
      _setstate(file_open_bit);   // Otherwise, file is open
      strcpy(_filepath, path);    // so save file name
   }
   set_error_on();
   return _MINCid;

}

//--------------------------------------------------------------------------
// MINC_File::close
// Close an open MINC volume. 
//--------------------------------------------------------------------------

void MINC_File::close(void) {
 
   set_error_off();
 
   if (this->is_open()){
      if (miclose(_MINCid) == MI_ERROR){
         _set_failure();            // could not close the file
      } else {
         _clrstate(file_open_bit);  // set status to closed file
         strcpy(_filepath, "");     // erase path name
         _MINCid = MI_ERROR;        // and reset the MINC id
      }
   }

   set_error_on();
}

//--------------------------------------------------------------------------
// MINC_File::attach_icv
// Attach the icv to the MINC file, after all icv properties have been
// modified.
//--------------------------------------------------------------------------

int MINC_File::attach_icv(void) {
   int status;

   set_error_off();
   if ((status = miicv_attach(_icvid, _MINCid, 
                              ncvarid(_MINCid, MIimage))) == MI_ERROR){
      _set_failure();               // could not attach icv
   } else {
      _setstate(icv_attached_bit);  // icv is attached
   }
   set_error_on();

   return status;
}

//--------------------------------------------------------------------------
// MINC_File::detach_icv
// Detachs an icv from the MINC volume.
//--------------------------------------------------------------------------

int MINC_File::detach_icv(void){
   _clrstate(icv_attached_bit);
   return miicv_detach(_icvid);
}

//--------------------------------------------------------------------------
// MINC_File::get_volume_info
// Returns information about the MINC volume.
//--------------------------------------------------------------------------

void MINC_File::get_volume_info(Volume_Info& volume_info) const {

   int i;

#ifdef DEBUG
   assert(this->is_good());
#endif

   volume_info.number_of_dimensions = _volume_info.number_of_dimensions;
   volume_info.datatype             = _volume_info.datatype;
   volume_info.valid_range[0]       = _volume_info.valid_range[0];
   volume_info.valid_range[1]       = _volume_info.valid_range[1];

   strcpy(volume_info.signtype, _volume_info.signtype);

   for (i=0; i < volume_info.number_of_dimensions; i++){
      volume_info.length[i] = _volume_info.length[i];
      volume_info.step[i]   = _volume_info.step[i];
      volume_info.start[i]  = _volume_info.start[i];
      strcpy(volume_info.dimension_names[i], _volume_info.dimension_names[i]);
   }
}

//--------------------------------------------------------------------------
// MINC_File::has_same_dimension_info_as
// Returns TRUE if minc_file has the same start/step information.
//--------------------------------------------------------------------------

int MINC_File::has_same_dimension_info_as(const MINC_File& minc_file) const {

   int same_info = TRUE;
   Volume_Info volume_info;

#ifdef DEBUG
   assert(this->is_good());
   assert(minc_file.is_good());
#endif
  
   minc_file.get_volume_info(volume_info);

   // Check if volumes have the same number of dimensions
   if (_volume_info.number_of_dimensions != volume_info.number_of_dimensions)
      same_info = FALSE;

   // Assume number of dimensions are the same.
   // Even if they are not previous test will have failed. 
   int n;
   for (n=0; n<_volume_info.number_of_dimensions; n++){ 
      if (_volume_info.step[n] != volume_info.step[n])
         same_info = FALSE;
      if (_volume_info.start[n] != volume_info.start[n])
         same_info = FALSE;
   }

   return same_info;
}

//--------------------------------------------------------------------------
// MINC_File::has_same_dimension_names_as
// Returns TRUE if minc_file has the same dimension names.
//--------------------------------------------------------------------------

int MINC_File::has_same_dimension_names_as(const MINC_File& minc_file) const {

   int same_names = TRUE;
   Volume_Info volume_info;

#ifdef DEBUG
   assert(this->is_good());
   assert(minc_file.is_good());
#endif
 
   minc_file.get_volume_info(volume_info);

   // Check if volumes have the same number of dimensions
   if (_volume_info.number_of_dimensions != volume_info.number_of_dimensions)
      same_names = FALSE;

   // Assume number of dimensions are the same.
   // Even if they are not previous test will have failed.
   int n;
   for (n=0; n<_volume_info.number_of_dimensions; n++){
      if (strcmp(_volume_info.dimension_names[n],
                 volume_info.dimension_names[n]) != 0)
         same_names = FALSE;
   }

   return same_names;

}

//--------------------------------------------------------------------------
// MINC_File::display_volume_info
// Writes formatted information about the MINC volume to an output stream.
//--------------------------------------------------------------------------

void MINC_File::display_volume_info(ostream& stream) const {

#ifdef DEBUG
   assert(this->is_good());
#endif

   // Output file information
   stream << "file: " << _filepath << endl;
   stream << "image: " << _volume_info.signtype << " "
          << MINC_File::type[_volume_info.datatype] << " "
          << _volume_info.valid_range[0] << " to "
          << _volume_info.valid_range[1] << endl;

   // Output dimension information
   stream << "image dimensions: ";
   int i;
   for(i=0; i< _volume_info.number_of_dimensions; i++){
      stream << _volume_info.dimension_names[i] << " ";
   }
   stream << endl;

   stream << setw(20) << "dimension name" << setw(8) << "length"
          << setw(12) << "step" << setw(12) << "start" << endl;
   stream << setw(20) << "--------------" << setw(8) << "------"
          << setw(12) << "----" << setw(12) << "-----" << endl;

   for(i=0; i< _volume_info.number_of_dimensions; i++){
      stream << setw(20) << _volume_info.dimension_names[i]
             << setw(8)  << _volume_info.length[i];

      if (_volume_info.step[i] == HUGE_VAL)
         stream << setw(12) << "unknown";
      else
         stream << setw(12) << _volume_info.step[i];

      if (_volume_info.start[i] == HUGE_VAL)
         stream << setw(12) << "unknown" << endl;
      else
         stream << setw(12) << _volume_info.start[i] << endl;

   }
   stream << endl;

}

//--------------------------------------------------------------------------
// MINC_File::put_attribute
// Sets an attribute in the MINC volume.  Overloaded for the required
// data type.
//--------------------------------------------------------------------------

int MINC_File::put_attribute(const char *varname, 
                             const char *attname, int value) {
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return miattputint(_MINCid, ncvarid(_MINCid, (char *)varname), 
                      (char *)attname, value);
}
      
int MINC_File::put_attribute(const char *varname, 
                             const char *attname, double value) {
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return miattputdbl(_MINCid, ncvarid(_MINCid, (char *)varname), 
                      (char *)attname, value);
}

int MINC_File::put_attribute(const char *varname, 
                             const char *attname, char *value) {
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return miattputstr(_MINCid, ncvarid(_MINCid, (char *)varname), 
                      (char *)attname, value);
}

//--------------------------------------------------------------------------
// MINC_File::get_attribute 
// Gets the value for a specified MINC attribute.  Overloaded for the
// required data type.
//--------------------------------------------------------------------------

int MINC_File::get_attribute(const char *varname, const char *attname, 
                             nc_type datatype,
                             int max_length, void *value, 
                             int *att_length) const {
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return miattget(_MINCid, ncvarid(_MINCid, (char *)varname), 
                   (char *)attname, datatype, max_length, 
                   value, att_length);
}

int MINC_File::get_attribute(const char *varname, const char *attname, 
                             nc_type datatype, void *value) const {
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return miattget1(_MINCid, ncvarid(_MINCid, (char *)varname), 
                    (char *)attname, datatype, value);
}

char *MINC_File::get_attribute(const char *varname, const char *attname, 
                               int maxlen, char *value) const {
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return miattgetstr(_MINCid, ncvarid(_MINCid, (char *)varname), 
                      (char *)attname, maxlen, value);
}

//--------------------------------------------------------------------------
// MINC_File::put_variable
// Sets the named MINC variable.  Overloaded for vector and scalar values.
//--------------------------------------------------------------------------

int MINC_File::put_variable(const char *varname, 
                            long start[], long count[],
                            nc_type datatype, char *sign, void *values){
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return mivarput(_MINCid, ncvarid(_MINCid, (char *)varname), 
                   start, count, datatype, sign, values);
}

int MINC_File::put_variable(const char *varname, 
                            long mindex[], nc_type datatype,
                            char *sign, void *value) {
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return mivarput1(_MINCid, ncvarid(_MINCid, (char *)varname), 
                    mindex, datatype, sign, value);
}

//--------------------------------------------------------------------------
// MINC_File::get_variable_type
// Returns the data type of the MINC volume.  Overloaded for enum
// and string values.
//--------------------------------------------------------------------------

int MINC_File::get_variable_type(const char *varname, char *type_name) const {
  
   int status;
   nc_type datatype;

   status = ncvarinq(_MINCid, ncvarid(_MINCid, (char *)varname), NULL, 
                  &datatype, NULL, NULL, NULL);

   switch(datatype){
      case NC_BYTE:
         strcpy(type_name,"byte");
         break;
      case NC_CHAR:
         strcpy(type_name,"char");
         break;
      case NC_SHORT:
         strcpy(type_name,"short");
         break;
      case NC_LONG:
         strcpy(type_name,"long");
         break;
      case NC_FLOAT:
         strcpy(type_name,"float");
         break;
      case NC_DOUBLE:
         strcpy(type_name,"double");
         break;
#if 0                           // These are obsolete.
      case NC_BITFIELD:
      case NC_STRING:
      case NC_IARRAY:
      case NC_DIMENSION:
      case NC_VARIABLE:
      case NC_ATTRIBUTE:
      case NC_UNSPECIFIED:
#endif
      default:
         strcpy(type_name,"unspecified or private type");
         break;
   }

   return status;
}

//--------------------------------------------------------------------------
// MINC_File::get_variable
// Gets the value for the specified MINC variable.  Overloaded for vector
// and scalar values. 
//--------------------------------------------------------------------------

int MINC_File::get_variable(const char *varname, 
                            long start[], long count[],
                            nc_type datatype, char *sign, 
                            void *values) const {
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return mivarget(_MINCid, ncvarid(_MINCid, (char *)varname), 
                   start, count, datatype, sign, values);
}

int MINC_File::get_variable(const char *varname, 
                            long mindex[], nc_type datatype,
                            char *sign, void *value) const {
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return mivarget1(_MINCid, ncvarid(_MINCid, (char *)varname), 
                    mindex, datatype, sign, value);
}

//--------------------------------------------------------------------------
// MINC_File::create_std_variable
// Creates a standard MINC variable in the volume.
//--------------------------------------------------------------------------

int MINC_File::create_std_variable(const char *varname, nc_type datatype,
                              int ndims, int dim[]) {
#ifdef DEBUG
   assert(!this->is_bad());
#endif
   return micreate_std_variable(_MINCid, (char *)varname, 
                                datatype, ndims, dim);
}

