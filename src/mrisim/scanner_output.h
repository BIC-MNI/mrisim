#ifndef __SCANNER_OUTPUT_H
#define __SCANNER_OUTPUT_H

//==========================================================================
// SCANNER_OUTPUT.H
// Scanner_Output class.
// Inherits from:  
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/scanner_output.h,v 1.1 2003-05-30 16:43:12 bert Exp $
 * $Log: scanner_output.h,v $
 * Revision 1.1  2003-05-30 16:43:12  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  16:08:30  rkwan
 * Release 3.1 - initial release.
 *
 *========================================================================*/

#include <minc/mriimage.h>
#include <minc/omincfile.h>

#include "mrisimargs.h"
#include "mriscanner.h"

//--------------------------------------------------------------------------
// Scanner_Output class
//--------------------------------------------------------------------------

class Scanner_Output {
   public:
      Scanner_Output(const mrisimArgs &args, char *time_stamp, 
                     MRI_Scanner &scanner);

      virtual ~Scanner_Output();

      inline int is_good(void) const;

      void display_info(const mrisimArgs &args, 
                        char *time_stamp,
                        const MRI_Scanner &scanner) const;
      void save_images(const mrisimArgs &args, MRI_Scanner &scanner);

      enum Output_Type {IMAGE_R, IMAGE_I, IMAGE_M, IMAGE_P,
                        RAW_R,   RAW_I,   RAW_M,   RAW_P};

   protected:
      int file_exists(const char *path) const;
      int create_output_file(const char *path, int clobber,
                      const char *time_stamp, const MRI_Scanner &scanner,
                      O_MINC_File &minc_file);
      char *extend_path(const char *original_path, int extension_length,
                      const char *extend_string) const;


   private:
      
      // --- Reconstructed image information --- //

      MRI_Image   _image;
      O_MINC_File _output;
      Output_Type  _output_type;

      // --- Raw data information --- //

      O_MINC_File _real_raw_data;
      O_MINC_File _imag_raw_data;

      int _good;
};

//--------------------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// Scanner_Output::is_good
// Returns TRUE if the Scanner_Output object is usable, i.e. all output
// files have been successfully opened.
//--------------------------------------------------------------------------

inline
int Scanner_Output::is_good(void) const {
   return _good;
}

#endif

