#ifndef __IOMINCFILE_H
#define __IOMINCFILE_H

//==========================================================================
// IOMINCFILE.H
// MINC input and output file class IO_MINC_File.
// Inherits from:  I_MINC_File, O_MINC_File
// Base class to:
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/iomincfile.h,v 1.1 2003-05-30 16:43:09 bert Exp $
 * $Log: iomincfile.h,v $
 * Revision 1.1  2003-05-30 16:43:09  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:27:15  rkwan
 * Release 2.5
 *
 * Revision 1.3  1995/12/11  15:13:07  rkwan
 * Fix RCS header bug.
 *
 *========================================================================*/

#include "imincfile.h"
#include "omincfile.h"

//--------------------------------------------------------------------------
// IO_MINC_File class
// Allows input and output of a MINC volume.
//--------------------------------------------------------------------------

class IO_MINC_File : public I_MINC_File, public O_MINC_File {
   public:
      IO_MINC_File();

      virtual ~IO_MINC_File();

      int open(const char *path);
};

#endif

