#ifndef __SPOILED_FLASH_H
#define __SPOILED_FLASH_H

/*****************************************************************************
 * 
 * SPOILED_FLASH.H
 *
 * Spoiled_FLASH Class
 *
 * R. Kwan
 * February 27, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/spoiled_flash.h,v 1.1 2003-05-30 16:43:14 bert Exp $
 * $Log: spoiled_flash.h,v $
 * Revision 1.1  2003-05-30 16:43:14  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:33:53  rkwan
 * Release 2.5
 *
 * Revision 2.1  1995/08/08  14:51:45  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "ffe.h"

/*****************************************************************************
 * Spoiled_FLASH Class
 *****************************************************************************/

class Spoiled_FLASH : public FFE {
   public:
      Spoiled_FLASH() : FFE() {}
      Spoiled_FLASH(Time_ms TR, Time_ms TE, Degrees flip_angle) :
           FFE(TR, TE, flip_angle) {}
      Spoiled_FLASH(Time_ms TR, Time_ms TE, Degrees flip_angle,
                    void (*f)(Vector_3D&,void *), void *data) :
           FFE(TR, TE, flip_angle, f, data) {}

      virtual ~Spoiled_FLASH() {}

      virtual void apply(Quick_Model& model);
      virtual void display_info(ostream& stream) const;
};

#endif
