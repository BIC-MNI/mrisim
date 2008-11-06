#ifndef __IMINCFILE_H
#define __IMINCFILE_H

//==========================================================================
// IMINCFILE.H
// MINC input file class I_MINC_File.
// Inherits from:  MINC_File
// Base class to:  IO_MINC_File
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/imincfile.h,v 1.3 2008-11-06 10:58:22 rotor Exp $
 * $Log: imincfile.h,v $
 * Revision 1.3  2008-11-06 10:58:22  rotor
 *  * fixed includes for iostream and friends
 *  * updated for new release (1.0.2)
 *
 * Revision 1.2  2004/08/10 15:19:59  bert
 * Add 'class' keyword
 *
 * Revision 1.1  2003/05/30 16:43:08  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:26:40  rkwan
 * Release 2.5
 *
 * Revision 1.3  1995/12/11  15:11:55  rkwan
 * Fix RCS header bug.
 *
 *========================================================================*/

#include <iostream>
#include "mincfile.h"

class O_MINC_File;
class MRI_Image;

//--------------------------------------------------------------------------
// I_MINC_File class
// Allows input of a MINC volume from disk.
//--------------------------------------------------------------------------

class I_MINC_File : virtual public MINC_File {
   friend class O_MINC_File;
   public:
      I_MINC_File();

      virtual ~I_MINC_File();

      // File creation routines
      virtual int open(const char *path);
      int attach_icv(void);

      // Data reading routines
      int load_hyperslab(long start[], long count[], void *volume);
      int load_slice(int slice_num, void *slice);
      int load_slice(int slice_num, MRI_Image& image);

   protected:

      void _get_volume_info(void);

};

#endif   
