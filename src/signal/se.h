#ifndef __SE_H
#define __SE_H

/*****************************************************************************
 * 
 * SE.H
 *
 * SE Class
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/se.h,v 1.1 2003-05-30 16:43:14 bert Exp $
 * $Log: se.h,v $
 * Revision 1.1  2003-05-30 16:43:14  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:33:37  rkwan
 * Release 2.5
 *
 * Revision 2.1  1995/08/08  14:46:20  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "quickseq.h"

/*****************************************************************************
 * SE Class
 *****************************************************************************/

class SE : public Quick_Sequence {
   public:
      SE();                            
      SE(Time_ms TR, Time_ms TE);
      SE(Time_ms TR, Time_ms TE, void (*f)(Vector_3D&,void *), void *data);
 
      virtual ~SE();
      virtual void apply(Quick_Model& model);
      virtual void display_info(ostream& stream) const;
      
   protected:
      Time_ms _TR;        // Repetition time
      Time_ms _TE;        // Echo time

};

#endif
