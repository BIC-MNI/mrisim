#ifndef __FISP_H
#define __FISP_H

/*****************************************************************************
 * 
 * FISP.H
 *
 * FISP Class
 *
 * R. Kwan
 * February 27, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/fisp.h,v 1.1 2003-05-30 16:43:13 bert Exp $
 * $Log: fisp.h,v $
 * Revision 1.1  2003-05-30 16:43:13  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:32:09  rkwan
 * Release 2.5
 *
 * Revision 2.1  1995/08/08  14:35:08  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "ffe.h"

/*****************************************************************************
 * FISP Class
 *****************************************************************************/

class FISP : public FFE {
   public:
      FISP() : FFE() {}
      FISP(Time_ms TR, Time_ms TE, Degrees flip_angle) :
           FFE(TR, TE, flip_angle) {}
      FISP(Time_ms TR, Time_ms TE, Degrees flip_angle,
                    void (*f)(Vector_3D&,void *), void *data) :
           FFE(TR, TE, flip_angle, f, data) {}

      virtual ~FISP() {}

      virtual void apply(Quick_Model& model);
      virtual void display_info(ostream& stream) const;
};

#endif
