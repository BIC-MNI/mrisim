#ifndef __FFE_H
#define __FFE_H

/*****************************************************************************
 * 
 * FFE.H
 *
 * FFE Class
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/ffe.h,v 1.1 2003-05-30 16:43:13 bert Exp $
 * $Log: ffe.h,v $
 * Revision 1.1  2003-05-30 16:43:13  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.1  1995/08/08  14:34:15  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "quickseq.h"

/*****************************************************************************
 * FFE Class
 *****************************************************************************/

class FFE : public Quick_Sequence {
   public:
      FFE();                            
      FFE(Time_ms TR, Time_ms TE, Degrees flip_angle);
      FFE(Time_ms TR, Time_ms TE, Degrees flip_angle,
         void (*f)(Vector_3D&,void *), void *data);
 
      virtual ~FFE();

      virtual void apply(Quick_Model& model) = 0;
      
   protected:
      Time_ms _TR;          // Repetition time
      Time_ms _TE;          // Echo time
      Degrees _flip_angle;

};

#endif
