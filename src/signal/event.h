#ifndef __EVENT_H
#define __EVENT_H

/*****************************************************************************
 * 
 * EVENT.H
 *
 * Event Class -- Base class for pulse sequence events.
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/event.h,v 1.2 2004-08-10 15:39:41 bert Exp $
 * $Log: event.h,v $
 * Revision 1.2  2004-08-10 15:39:41  bert
 * Add 'class' keyword to friend declaration
 *
 * Revision 1.1  2003/05/30 16:43:12  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.2  1995/12/11  14:23:12  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.1  1995/08/08  14:33:23  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include <minc/mristring.h>
#include "vector.h"
#include "spin_model.h"

// Event types
#define EMPTY_EVENT 0
#define RF_PULSE    1
#define SAMPLE      2
#define SPOILER     3
#define REPEAT      4

/*****************************************************************************
 * Event Class
 *****************************************************************************/
class Custom_Sequence;

class Event {
   public:
      Event(); 
      Event(Time_ms t);
      Event(const Event& event);

      virtual ~Event() {}

      virtual Vector_3D& apply(Spin_Model& m) = 0;        // apply event to mag.
      virtual void get_descriptor_string(char s[]);
      virtual Event *make_new_copy_of_event(void) = 0;
      
      Time_ms get_event_time(void);
      int     get_event_id(void);
      int     get_event_type(void);

      void    set_event_time(Time_ms t);
      void    set_event_id(int event_id);
      void    set_event_type(int event_type);

      friend  class Custom_Sequence;

   protected:
      Time_ms _event_time;                 // time stamp of event
      int     _event_id;                   // unique event identifier
      int     _event_type;                 // type of event

      Event   *_next;                      // next event in linked list
      static int _last_id;
};

#endif

