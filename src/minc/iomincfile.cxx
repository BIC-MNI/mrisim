//==========================================================================
// IOMINCFILE.CXX
// MINC input and output file class IO_MINC_File.
// Inherits from:  I_MINC_File, O_MINC_File
// Base class to:
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/iomincfile.cxx,v 1.1 2003-05-30 16:43:08 bert Exp $
// $Log: iomincfile.cxx,v $
// Revision 1.1  2003-05-30 16:43:08  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:27:22  rkwan
// Release 2.5
//
// Revision 1.2  1995/12/11  14:15:59  rkwan
// Updated for fast_iso_model.
//
//==========================================================================

#include "iomincfile.h"
#include <string.h>

//--------------------------------------------------------------------------
// IO_MINC_File constructor
//--------------------------------------------------------------------------

IO_MINC_File::IO_MINC_File() : I_MINC_File(), O_MINC_File() {}

//--------------------------------------------------------------------------
// IO_MINC_File destructor
//--------------------------------------------------------------------------

IO_MINC_File::~IO_MINC_File() {}

//--------------------------------------------------------------------------
// IO_MINC_File::open
// Opens a MINC file for reading and writing.
// Sets status bits on failure, and returns the id of the open MINC file.
//--------------------------------------------------------------------------

int IO_MINC_File::open(const char *path) {

   // If file is already open close it
   if (this->is_open()){
      this->close();
   }
   
   // Open the file
   set_error_off();
   if ((_MINCid = miopen((char *)path, NC_WRITE)) == MI_ERROR){
      _set_failure();              // Could not open file
   } else {
      _setstate(file_open_bit);    // Otherwise, file is open
      strcpy(_filepath, path);     // so save file name
   }
   set_error_on();

   return _MINCid;

}
