#ifndef __FLASH_H
#define __FLASH_H

/*****************************************************************************
 * 
 * FLASH.H
 *
 * FLASH Class
 *
 * R. Kwan
 * February 27, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/flash.h,v 1.1 2003-05-30 16:43:13 bert Exp $
 * $Log: flash.h,v $
 * Revision 1.1  2003-05-30 16:43:13  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:32:30  rkwan
 * Release 2.5
 *
 * Revision 2.1  1995/08/08  14:35:52  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "ffe.h"

/*****************************************************************************
 * FLASH Class
 *****************************************************************************/

class FLASH : public FFE {
   public:
      FLASH() : FFE() {}
      FLASH(Time_ms TR, Time_ms TE, Degrees flip_angle) :
           FFE(TR, TE, flip_angle) {}
      FLASH(Time_ms TR, Time_ms TE, Degrees flip_angle,
                    void (*f)(Vector_3D&,void *), void *data) :
           FFE(TR, TE, flip_angle, f, data) {}

      virtual ~FLASH() {}

      virtual void apply(Quick_Model& model);
      virtual void display_info(ostream& stream) const;
};

#endif
