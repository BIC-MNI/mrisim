#ifndef __REPEAT_H
#define __REPEAT_H

/*****************************************************************************
 * 
 * REPEAT.H
 *
 * Repeat Class
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/repeat.h,v 1.1 2003-05-30 16:43:13 bert Exp $
 * $Log: repeat.h,v $
 * Revision 1.1  2003-05-30 16:43:13  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.2  1995/12/11  14:26:54  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.1  1995/08/08  14:44:05  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "event.h"

/*****************************************************************************
 * Repeat Class
 *****************************************************************************/

class Repeat : public Event {
   public:
      Repeat() : Event() {}            // constructors
      Repeat(Time_ms t) : Event(t) {}
      Repeat(const Repeat& repeat) : Event((Event &)repeat) {}

      virtual Vector_3D& apply(Spin_Model& m);    
      virtual void get_descriptor_string(char s[]);
      virtual Event *make_new_copy_of_event(void);
  
};

#endif

