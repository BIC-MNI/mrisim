#ifndef __IR_H
#define __IR_H

/*****************************************************************************
 * 
 * IR.H
 *
 * IR Class
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/ir.h,v 1.1 2003-05-30 16:43:13 bert Exp $
 * $Log: ir.h,v $
 * Revision 1.1  2003-05-30 16:43:13  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:32:44  rkwan
 * Release 2.5
 *
 * Revision 2.1  1995/08/08  14:36:39  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "quickseq.h"

/*****************************************************************************
 * IR Class
 *****************************************************************************/

class IR : public Quick_Sequence {
   public:
      IR();                            
      IR(Time_ms TR, Time_ms TE, Time_ms TI);
      IR(Time_ms TR, Time_ms TE, Time_ms TI,
         void (*f)(Vector_3D&,void *), void *data);
 
      virtual ~IR();
      virtual void apply(Quick_Model& mag);
      virtual void display_info(ostream& stream) const;
      
   protected:
      Time_ms _TR;        // Repetition time
      Time_ms _TE;        // Echo time
      Time_ms _TI;        // Inversion time

};

#endif
