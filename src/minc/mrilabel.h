#ifndef __MRILABEL_H
#define __MRILABEL_H

//===========================================================================
// MRILABEL.H
// Labelled slice class.
// Inherits from:  MRI_Byte_Matrix
// Base class to:  
// 
// R.Kwan
// April, 19, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

/*===========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/mrilabel.h,v 1.1 2003-05-30 16:43:09 bert Exp $
 * $Log: mrilabel.h,v $
 * Revision 1.1  2003-05-30 16:43:09  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:48:05  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:23:02  rkwan
 * Release 2.5
 *
 *
 *=========================================================================*/

#include <limits.h>
#include "mrimatrix.h"

typedef unsigned char Label;

#ifndef MAX_LABEL
#define MAX_LABEL UCHAR_MAX
#endif

//---------------------------------------------------------------------------
// MRI_Label class
// Specialization of MRI_Byte_Matrix to store slices of a labelled volume.
//---------------------------------------------------------------------------

class MRI_Label : public MRI_Byte_Matrix {
   public:
      MRI_Label(unsigned int nrows, unsigned int ncols) :
         MRI_Byte_Matrix(nrows, ncols) {}
      MRI_Label(const MRI_Byte_Matrix& mat) : MRI_Byte_Matrix(mat) {}
      MRI_Label(const MRI_Short_Matrix& mat) : MRI_Byte_Matrix(mat) {}
      MRI_Label(const MRI_Float_Matrix& mat) : MRI_Byte_Matrix(mat) {}
      MRI_Label(const MRI_Double_Matrix& mat) : MRI_Byte_Matrix(mat) {}

      virtual ~MRI_Label() {}

      void get_mask(Label label, MRI_Byte_Matrix& mask) const;
      void get_mask(Label label, MRI_Short_Matrix& mask) const;
      void get_mask(Label label, MRI_Float_Matrix& mask) const;
      void get_mask(Label label, MRI_Double_Matrix& mask) const;

   private:
      MRI_Label() : MRI_Byte_Matrix() {}
};

#endif

