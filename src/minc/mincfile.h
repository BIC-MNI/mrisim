#ifndef __MINCFILE_H
#define __MINCFILE_H

//==========================================================================
// MINCFILE.H
// Base class for MINC file classes.
// Inherits from:  MINC_ICV
// Base class to:  I_MINC_File, O_MINC_File
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/mincfile.h,v 1.2 2008-11-06 10:58:23 rotor Exp $
 * $Log: mincfile.h,v $
 * Revision 1.2  2008-11-06 10:58:23  rotor
 *  * fixed includes for iostream and friends
 *  * updated for new release (1.0.2)
 *
 * Revision 1.1  2003/05/30 16:43:09  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:21:59  rkwan
 * Release 2.5
 *
 * Revision 1.4  1996/01/17  17:59:21  rkwan
 * Update for tx/rx coil modelling.
 *
 * Revision 1.3  1995/12/11  15:11:08  rkwan
 * Fix RCS header bug.
 *
 *========================================================================*/

extern "C" {
#include <minc.h>
#undef public
}
#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include "mincicv.h"
#include "mrimatrix.h"
#include "mrivolume.h"

//--------------------------------------------------------------------------
// Volume dimensions
//--------------------------------------------------------------------------

#ifndef SLICE
#define SLICE 0
#define ROW   1
#define COLUMN 2
#endif

//--------------------------------------------------------------------------
// Volume information structure
//--------------------------------------------------------------------------

typedef struct {
   int number_of_dimensions;                         // num of dimensions
   int axes[MAX_VAR_DIMS];                           // axis tranformation
   long length[MAX_VAR_DIMS];                        // num of elements
   double step[MAX_VAR_DIMS];                        // step size
   double start[MAX_VAR_DIMS];                       // start of dimension
   char dimension_names[MAX_VAR_DIMS][MAX_NC_NAME];  // name of dimension
   nc_type datatype;                                 // volume data type
   char signtype[MI_MAX_ATTSTR_LEN];                 // signed or unsigned?
   double valid_range[2];                            // valid data range
} Volume_Info;
   
//--------------------------------------------------------------------------
// MINC_File class
// Base class for MINC file classes.
//--------------------------------------------------------------------------

class MINC_File : public MINC_ICV {
   public:
      MINC_File();
      virtual ~MINC_File();

      // --- General file operations --- //

      virtual int open(const char *path, int mode);
      virtual int attach_icv(void);
      void        close(void);
      int         detach_icv(void);

      // --- Query object state --- //

      inline int is_good(void) const;
      inline int is_bad(void) const;
      inline int is_open(void) const;
      inline int has_ICV_attached(void) const;

      // --- Error Handling --- //

      inline void set_error_on(void);
      inline void set_error_off(void);
      inline void set_full_error(void);
      inline void set_verbose_error(void);
      inline void set_fatal_error(void);
      inline int  get_error_status(void) const;

      // --- Volume information convenience functions --- //
    
      inline const char *get_filename(void) const;
      inline int         get_nrows(void) const;
      inline int         get_ncols(void) const;
      inline int         get_nslices(void) const;
      inline int         get_num_of_slice_elements(void) const;
      inline double      get_step(int n) const;
      inline double      get_start(int n) const;

      void get_volume_info(Volume_Info &volume_info) const;
      void display_volume_info(ostream& stream) const;

      inline int is_same_slice_size_as(const MRI_Matrix& mat) const;
      inline int is_same_size_as(const MRI_Volume& vol) const;
      inline int is_same_size_as(const MINC_File& minc_file) const;

      int has_same_dimension_info_as(const MINC_File& minc_file) const;
      int has_same_dimension_names_as(const MINC_File& minc_file) const;
      
      // --- Low-level NetCDF and MINC functions --- //

      inline int get_cdfid(void) const;
      inline int set_fill_on(void); 
      inline int set_fill_off(void);
      inline int abort(void);
      inline int sync(void);
  
      // --- MINC attribute access functions --- //

      int put_attribute(const char *varname, 
                        const char *attname, int value);
      int put_attribute(const char *varname, 
                        const char *attname, double value);
      int put_attribute(const char *varname, 
                        const char *attname, char *value);

      int get_attribute(const char *varname, const char *attname, 
                        nc_type datatype, int max_length, 
                        void *value, int *att_length) const;
      int get_attribute(const char *varname, const char *attname, 
                        nc_type datatype, void *value) const;
      char *get_attribute(const char *varname, const char *attname, 
                        int maxlen, char *value) const;

      // --- MINC variable access functions --- //

      int put_variable(const char *varname, long start[], 
                       long count[], nc_type datatype, 
                       char *sign, void *values);
      int put_variable(const char *varname, long mindex[], 
                       nc_type datatype,
                       char *sign, void *value);

      inline int     get_variable_id(const char *varname) const;
      inline nc_type get_variable_type(const char *varname) const;

      int get_variable_type(const char *varname, char *type_name) const;

      int get_variable(const char *varname, long start[], 
                       long count[], nc_type datatype, 
                       char *sign, void *values) const;
      int get_variable(const char *varname, long mindex[], 
                       nc_type datatype, char *sign, void *value) const;

      int create_std_variable(const char *varname, nc_type datatype,
                       int ndims, int dim[]);

   protected:

      // --- Internal data structures --- //

      int         _MINCid;        // MINC id number
      Volume_Info _volume_info;   // Volume information for the MINC volume
      char        _filepath[255]; // Path for the MINC volume on disk

      // Static look up tables
      static char *type[];
      static double default_min[][2];
      static double default_max[][2];

      // --- Error Handling --- //
      enum io_state {goodbit = 0, failbit = 1, 
                     file_open_bit = 2, icv_attached_bit = 4};

      void _setstate(int state) {_state |= state;}
      void _set_failure(void)   {_state |= failbit;}
      void _clrstate(int state) {_state &= ~state;}

   private:
 
      // --- Internal data structures --- //
      static int _ncoldopts;
      int        _state;             // Error state of the object

};

//---------------------------------------------------------------------------
// Inline member functions
//---------------------------------------------------------------------------

//--------------------------------------------------------------------------
// MINC_File::is_good
// Returns TRUE if a MINC volume is usable, i.e. a file is open and there
// is an attached icv.
//--------------------------------------------------------------------------

inline
int MINC_File::is_good(void) const {
   return (_state == (goodbit | file_open_bit | icv_attached_bit));
}

//--------------------------------------------------------------------------
// MINC_File::is_bad
// Returns TRUE if a MINC has failed on the MINC volume.
//--------------------------------------------------------------------------

inline
int MINC_File::is_bad(void) const {
   return (_state & failbit) || (_MINCid == MI_ERROR);
}

//--------------------------------------------------------------------------
// MINC_File::is_open
// Returns TRUE if the MINC volume is open.
//--------------------------------------------------------------------------

inline
int MINC_File::is_open(void) const {
   return _state & file_open_bit;
}

//--------------------------------------------------------------------------
// MINC_File::has_ICV_attached
// Returns TRUE if the MINC volume has an attached ICV associated
// with it.
//--------------------------------------------------------------------------

inline
int MINC_File::has_ICV_attached(void) const {
   return _state & icv_attached_bit;
}

//--------------------------------------------------------------------------
// MINC_File::set_error_on
// Turns the previous error checking options on.
//--------------------------------------------------------------------------

inline
void MINC_File::set_error_on(void) {
   set_ncopts(_ncoldopts);
}

//--------------------------------------------------------------------------
// MINC_File::set_error_off
// Turns error checking on.
//--------------------------------------------------------------------------

inline
void MINC_File::set_error_off(void) {
   _ncoldopts =get_ncopts();
   set_ncopts(0);
}

//--------------------------------------------------------------------------
// MINC_File::set_full_error
// Turns on full verbose and fatal error checking.
//--------------------------------------------------------------------------

inline
void MINC_File::set_full_error(void) {
   _ncoldopts =get_ncopts();
   set_ncopts(NC_VERBOSE | NC_FATAL);
}

//--------------------------------------------------------------------------
// MINC_File::set_verbose_error
// Turns on verbose error checking.
//--------------------------------------------------------------------------

inline 
void MINC_File::set_verbose_error(void) {
   _ncoldopts =get_ncopts();
   set_ncopts(ncopts | NC_VERBOSE);
}

//--------------------------------------------------------------------------
// MINC_File::set_fatal_error
// Turns on fatal error checking.
//--------------------------------------------------------------------------

inline
void MINC_File::set_fatal_error(void) {
   _ncoldopts =get_ncopts();
   set_ncopts(ncopts | NC_FATAL);
}

//--------------------------------------------------------------------------
// MINC_File::get_error_status
// Returns the current error status of the MINC library.
//--------------------------------------------------------------------------

inline
int  MINC_File::get_error_status(void) const {
   return ncerr;
}

//--------------------------------------------------------------------------
// MINC_File::get_filename
// Returns the stored file name of the MINC file.
//--------------------------------------------------------------------------

inline
const char *MINC_File::get_filename(void) const {
   return (const char *)_filepath;
}

//--------------------------------------------------------------------------
// MINC_File::get_nrows
// Returns the number of rows in the volume.
//--------------------------------------------------------------------------

inline
int MINC_File::get_nrows(void) const {
   return _volume_info.length[ROW];
}

//--------------------------------------------------------------------------
// MINC_File::get_ncols
// Returns the number of columns in the volume.
//--------------------------------------------------------------------------

inline
int MINC_File::get_ncols(void) const {
   return _volume_info.length[COLUMN];
}

//--------------------------------------------------------------------------
// MINC_File::get_nslices
// Returns the number of slices in the volume.
//--------------------------------------------------------------------------

inline
int MINC_File::get_nslices(void) const {
   return _volume_info.length[SLICE];
}

//--------------------------------------------------------------------------
// MINC_File::get_num_of_slice_elements
// Returns the number of elements in a slice.
//--------------------------------------------------------------------------

inline
int MINC_File::get_num_of_slice_elements(void) const {
   return _volume_info.length[ROW]*_volume_info.length[COLUMN];
}

//--------------------------------------------------------------------------
// MINC_File::get_step
// Returns the step size for the specified index.
//--------------------------------------------------------------------------

inline
double MINC_File::get_step(int n) const {
   return _volume_info.step[n]; 
}

//--------------------------------------------------------------------------
// MINC_File::get_start
// Returns the start coordinate for the specified index.
//--------------------------------------------------------------------------

inline
double MINC_File::get_start(int n) const {
   return _volume_info.start[n]; 
}

//--------------------------------------------------------------------------
// MINC_File::is_same_slice_size_as
// Returns TRUE if the matrix has the same row and column dimensions as
// the MINC file.
//---------------------------------------------------------------------------

inline
int MINC_File::is_same_slice_size_as(const MRI_Matrix& mat) const {
   int same_size;

   if ((this->get_nrows() == mat.get_nrows()) &&
       (this->get_ncols() == mat.get_ncols())){
      same_size = TRUE;
   } else {
      same_size = FALSE;
   }
  
   return same_size;
}

//---------------------------------------------------------------------------
// MINC_File::is_same_size_as
// Returns TRUE if the volume has the same row, column, and slice dimensions
// as the MINC file.
//---------------------------------------------------------------------------

inline
int MINC_File::is_same_size_as(const MRI_Volume& vol) const {
   int same_size;

   if ((this->get_nrows() == vol.get_nrows()) &&
       (this->get_ncols() == vol.get_ncols()) &&
       (this->get_nslices() == vol.get_nslices())){
      same_size = TRUE;
   } else {
      same_size = FALSE;
   }
  
   return same_size;
}

//--------------------------------------------------------------------------
// MINC_File::is_same_size_as
// Returns TRUE if minc_file has the same dimension size information.
//--------------------------------------------------------------------------

inline
int MINC_File::is_same_size_as(const MINC_File& minc_file) const {

   int same_size = TRUE;

   if ((get_nrows() == minc_file.get_nrows()) &&
       (get_ncols() == minc_file.get_ncols()) &&
       (get_nslices() == minc_file.get_nslices())){
      same_size = TRUE;
   } else {
      same_size = FALSE;
   }
   return same_size;

}

//---------------------------------------------------------------------------
// MINC_File::get_cdfid
// Return the id of the MINC_File.
//---------------------------------------------------------------------------

inline
int MINC_File::get_cdfid(void) const {
   return _MINCid;
}

//---------------------------------------------------------------------------
// MINC_File::set_fill_on
// Turns on variable filling. 
//---------------------------------------------------------------------------

inline
int MINC_File::set_fill_on(void) {
   return ncsetfill(_MINCid, NC_FILL);
}

//---------------------------------------------------------------------------
// MINC_File::set_fill_off
//---------------------------------------------------------------------------

inline
int MINC_File::set_fill_off(void) {
   return ncsetfill(_MINCid, NC_NOFILL);
}

//---------------------------------------------------------------------------
// MINC_File::abort
// Abort any changes to the file if it is in define mode.
//---------------------------------------------------------------------------

inline
int MINC_File::abort(void) {
   return ncabort(_MINCid);
}

//---------------------------------------------------------------------------
// MINC_File::sync
// Synchronize an open file to disk.
//---------------------------------------------------------------------------

inline
int MINC_File::sync(void) {
   return ncsync(_MINCid);
}

//---------------------------------------------------------------------------
// MINC_File::get_variable_id
// Returns the id of the named MINC variable.
//---------------------------------------------------------------------------

inline
int MINC_File::get_variable_id(const char *varname) const {
   return ncvarid(_MINCid, (char *)varname);
} 

//--------------------------------------------------------------------------
// MINC_File::get_variable_type
// Returns the data type of the MINC volume.  Overloaded for enum
// and string values.
//--------------------------------------------------------------------------

inline
nc_type MINC_File::get_variable_type(const char *varname) const {
   nc_type datatype;
   ncvarinq(_MINCid, ncvarid(_MINCid, (char *)varname), NULL,
            &datatype, NULL, NULL, NULL);
   return datatype;
}

#endif    
