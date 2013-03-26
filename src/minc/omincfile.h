#ifndef __OMINCFILE_H
#define __OMINCFILE_H

//==========================================================================
// OMINCFILE.H
// MINC output file class O_MINC_File.
// Inherits from:  MINC_File
// Base class to:  IO_MINC_File
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/omincfile.h,v 1.2 2004-08-10 15:35:55 bert Exp $
 * $Log: omincfile.h,v $
 * Revision 1.2  2004-08-10 15:35:55  bert
 * Add 'class' keyword to friend declarations
 *
 * Revision 1.1  2003/05/30 16:43:10  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:29:25  rkwan
 * Release 2.5
 *
 * Revision 1.4  1996/01/17  17:59:39  rkwan
 * Update for tx/rx coil modelling.
 *
 * Revision 1.3  1995/12/11  15:12:25  rkwan
 * Fix RCS header bug.
 *
 *========================================================================*/

#include "mincfile.h"
#include "mrimatrix.h"
#include "mriimage.h"

class I_MINC_File;

//--------------------------------------------------------------------------
// O_MINC_File class
// Allows output of a MINC volume to disk.
//--------------------------------------------------------------------------

class O_MINC_File : virtual public MINC_File {
   friend class I_MINC_File;
   public:
      O_MINC_File();

      virtual ~O_MINC_File();

      // File creation routines
      virtual int open(const char *path);
      int create(const char *path, int mode = NC_NOCLOBBER); 

      // Data writing routines
      int save_hyperslab(long start[], long count[], void *volume);
      int save_slice(int slice_num, void *slice);
      int save_slice(int slice_num, MRI_Image& image);

      // VIO_Volume information routines
      void set_volume_info(const I_MINC_File& ifile, 
                           const char *argstring = NULL);
      void set_volume_info(int cdfid, const Volume_Info &volume_info,
                           const char *argstring = NULL);

   protected:

      void _setup_image_variables(int inMINCid, int ndims, int dim[]);
      void _update_history(const char *arg_string);

};

#endif
