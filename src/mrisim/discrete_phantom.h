#ifndef __DISCRETE_PHANTOM_H
#define __DISCRETE_PHANTOM_H

//==========================================================================
// DISCRETE_PHANTOM.H
// Discrete_Phantom class.
// Inherits from:  Discrete_Label_Phantom, Tissue_Phantom
// Base class to:  
//
// R. Kwan
// (C) Copyright 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/discrete_phantom.h,v 1.1 2003-05-30 16:43:10 bert Exp $
 * $Log: discrete_phantom.h,v $
 * Revision 1.1  2003-05-30 16:43:10  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:09:07  rkwan
 * Release 2.5
 *
 *========================================================================*/

#include "discrete_label_phantom.h"
#include "tissue_phantom.h"

//---------------------------------------------------------------------------
// Discrete_Phantom class
// Describes a discrete labelled volume MRI phantom.
//---------------------------------------------------------------------------

class Discrete_Phantom : public Discrete_Label_Phantom,
                         public Tissue_Phantom {
   public:
      Discrete_Phantom(unsigned int n_tissue_classes); 

      virtual ~Discrete_Phantom();

      // --- Image Formation --- // 
      void get_simulated_phantom_slice(int slice_num,
                                       Complex_Slice& sim_slice);
      void get_simulated_mag_phantom_slice(int slice_num,
                                       Real_Slice& sim_slice);
      void get_simulated_real_phantom_slice(int slice_num,
                                       Real_Slice& sim_slice);

};

#endif

