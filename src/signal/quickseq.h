#ifndef __QUICKSEQ_H
#define __QUICKSEQ_H

/*****************************************************************************
 * 
 * QUICKSEQ.H
 *
 * Quick_Sequence class.
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/quickseq.h,v 1.1 2003-05-30 16:43:13 bert Exp $
 * $Log: quickseq.h,v $
 * Revision 1.1  2003-05-30 16:43:13  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.1  1995/08/08  14:42:58  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include <mrisim/mrisim.h>
#include "pulseseq.h"
#include "quick_model.h"
#include "sample.h"

/*****************************************************************************
 * Quick_Sequence Class
 *****************************************************************************/

class Quick_Sequence : public Pulse_Sequence {
   public:
      Quick_Sequence();
      Quick_Sequence(void (*f)(Vector_3D&,void *), void *data);
      Quick_Sequence(const Quick_Sequence& p);

      virtual ~Quick_Sequence();

      virtual void apply(Quick_Model& model) = 0;

   protected:
      Sample *_sample;
};

#endif
