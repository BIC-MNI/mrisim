#ifndef __CE_FAST_H
#define __CE_FAST_H

/*****************************************************************************
 * 
 * CE_FAST.H
 *
 * CE_FAST Class
 *
 * R. Kwan
 * February 27, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/ce_fast.h,v 1.1 2003-05-30 16:43:12 bert Exp $
 * $Log: ce_fast.h,v $
 * Revision 1.1  2003-05-30 16:43:12  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:31:19  rkwan
 * Release 2.5
 *
 * Revision 2.2  1995/08/08  14:27:52  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "ffe.h"

/*****************************************************************************
 * CE_FAST Class
 *****************************************************************************/

class CE_FAST : public FFE {
   public:
      CE_FAST() : FFE() {}
      CE_FAST(Time_ms TR, Time_ms TE, Degrees flip_angle) :
           FFE(TR, TE, flip_angle) {}
      CE_FAST(Time_ms TR, Time_ms TE, Degrees flip_angle,
                    void (*f)(Vector_3D&,void *), void *data) :
           FFE(TR, TE, flip_angle, f, data) {}

      virtual ~CE_FAST() {}

      virtual void apply(Quick_Model& model);
      virtual void display_info(ostream& stream) const;
};

#endif
