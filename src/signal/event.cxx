/*****************************************************************************
 * 
 * EVENT.CXX
 *
 * Event Class -- Base class for pulse sequence events.
 *
 * R. Kwan
 * August 3, 1994
 *
 * (C) Copyright 1994 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/event.cxx,v 1.1 2003-05-30 16:43:12 bert Exp $
// $Log: event.cxx,v $
// Revision 1.1  2003-05-30 16:43:12  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.2  1995/12/11  14:23:19  rkwan
// Updated for fast_iso_model.
//
// Revision 2.1  1995/08/08  14:33:43  rkwan
// Updated for new spin model and pulse sequence.
//

#include "event.h"

/*****************************************************************************
 * Event Class
 *****************************************************************************/

int Event::_last_id = 0;

Event::Event() {
   _event_time = 0;
   _event_id   = _last_id++;
   _event_type = EMPTY_EVENT;
   _next       = (Event *)NULL;
}

Event::Event(Time_ms t){
   _event_time = t;
   _event_id   = _last_id++;
   _event_type = EMPTY_EVENT;
   _next       = (Event *)NULL;
}

Event::Event(const Event& event){
   _event_time = event._event_time;
   _event_id   = _last_id++;
   _event_type = event._event_type;
   _next       = (Event *)NULL;
}

void Event::get_descriptor_string(char s[]){
   strcpy((char *)s, "EMPTY EVENT");
}

Time_ms Event::get_event_time(void){
   return _event_time;
}

int Event::get_event_id(void){
   return _event_id;
}

int Event::get_event_type(void){
   return _event_type;
}

void Event::set_event_time(Time_ms event_time){
   _event_time = event_time;
}

void Event::set_event_id(int event_id){
   _event_id = event_id;
}

void Event::set_event_type(int event_type){
   _event_type = event_type;
}

