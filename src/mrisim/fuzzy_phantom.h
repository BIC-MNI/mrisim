#ifndef __FUZZY_PHANTOM_H
#define __FUZZY_PHANTOM_H

//==========================================================================
// FUZZY_PHANTOM.H
// Fuzzy_Phantom class.
// Inherits from:  Fuzzy_Label_Phantom, Tissue_Phantom
// Base class to:  
//
// R. Kwan
// (C) Copyright 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/fuzzy_phantom.h,v 1.1 2003-05-30 16:43:10 bert Exp $
 * $Log: fuzzy_phantom.h,v $
 * Revision 1.1  2003-05-30 16:43:10  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:51:30  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:10:35  rkwan
 * Release 2.5
 *
 *========================================================================*/

#include "fuzzy_label_phantom.h"
#include "tissue_phantom.h"

//---------------------------------------------------------------------------
// Fuzzy_Phantom class
// Describes a fuzzy labelled MRI phantom.
//---------------------------------------------------------------------------

class Fuzzy_Phantom : public Fuzzy_Label_Phantom,
                      public Tissue_Phantom {
   public:
      Fuzzy_Phantom(unsigned int n_tissue_classes);

      virtual ~Fuzzy_Phantom();

      // --- Image Formation --- // 
      void get_simulated_phantom_slice(int slice_num,
                                       Complex_Slice& sim_slice);
      void get_simulated_mag_phantom_slice(int slice_num,
                                       Real_Slice& sim_slice);
      void get_simulated_real_phantom_slice(int slice_num,
                                       Real_Slice& sim_slice);
};

#endif
