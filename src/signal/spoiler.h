#ifndef __SPOILER_H
#define __SPOILER_H

/*****************************************************************************
 * 
 * SPOILER.H
 *
 * Spoiler Class
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/spoiler.h,v 1.1 2003-05-30 16:43:14 bert Exp $
 * $Log: spoiler.h,v $
 * Revision 1.1  2003-05-30 16:43:14  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.2  1995/12/11  14:28:35  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.1  1995/08/08  14:54:35  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "event.h"

/*****************************************************************************
 * Spoiler Class
 *****************************************************************************/

class Spoiler : public Event {
   public:
      Spoiler() : Event() {}           // constructors
      Spoiler(Time_ms t) : Event(t) {}
      Spoiler(const Spoiler& spoiler) : Event((Event&)spoiler) {}

      virtual Vector_3D& apply(Spin_Model& m);
      virtual void get_descriptor_string(char s[]);
      virtual Event *make_new_copy_of_event(void);
};

#endif

