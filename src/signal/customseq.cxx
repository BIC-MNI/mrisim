/*****************************************************************************
 * 
 * CUSTOMSEQ.CXX
 *
 * Custom_Sequence Class
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/customseq.cxx,v 1.1 2003-05-30 16:43:12 bert Exp $
// $Log: customseq.cxx,v $
// Revision 1.1  2003-05-30 16:43:12  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:31:36  rkwan
// Release 2.5
//
// Revision 2.2  1995/12/11  14:23:03  rkwan
// Updated for fast_iso_model.
//
// Revision 2.1  1995/08/08  14:32:27  rkwan
// Updated for new spin model and pulse sequence.
//

#include "customseq.h"
#include <stdio.h>

/*****************************************************************************
 * Custom_Sequence Class
 *****************************************************************************/

/*****************************************************************************
 * Custom_Sequence constructors
 *****************************************************************************/

Custom_Sequence::Custom_Sequence() : Pulse_Sequence() {

   _num_events      = 0;
   _event_list      = (Event *)NULL;
   _next_event      = _event_list;

}

Custom_Sequence::Custom_Sequence(const Pulse_Sequence& p) :
   Pulse_Sequence(p) {

   _num_events      = 0;
   _event_list      = (Event *)NULL;
   _next_event      = _event_list;

}

Custom_Sequence::Custom_Sequence(const Custom_Sequence& p) :
   Pulse_Sequence((Pulse_Sequence&)p) {

   _delete_list(_event_list);
   Event *source = p._event_list;
   _num_events = p._num_events;

   if (_num_events != 0){
      _event_list = source->make_new_copy_of_event();
      Event *target = _event_list;
      source = source->_next;

      while(source != NULL){
         target->_next = source->make_new_copy_of_event();
         target = target->_next;
         source = source->_next;
      }
      _next_event = _event_list;

   } else {
      _event_list = (Event *)NULL;
      _next_event = _event_list;
      _num_events = 0;
   }

}

/*****************************************************************************
 * Custom_Sequence destructor
 *****************************************************************************/

Custom_Sequence::~Custom_Sequence(){
   _delete_list(_event_list);
}

/*****************************************************************************
 * Custom_Sequence::add_event
 * Add an event to the pulse sequence event_list.
 * The event is added to the list according to its event_time, so that
 * the list remains ordered.   _current_event is set to point to the
 * newly added event.
 *****************************************************************************/

void Custom_Sequence::add_event(Event *event){

   if (_num_events == 0){

      // Simply add the event to an empty event list
      _event_list = event;
      _next_event = _event_list;
      _num_events++;

   } else {

      Event *ptr = _find_event_after_time(_next_event->_event_time);

      if (ptr == NULL){

         if (_next_event == NULL){

            // Event time is before the current first event
            ptr = _event_list;
            _event_list = event;
            event->_next = ptr;
            _next_event = _event_list;

         } else {

            // Event time is after the last event
            _next_event->_next = event;
            _next_event = _next_event->_next;
         }

      } else {
         
         // Insert the new event
         _next_event->_next = event;
         event->_next = ptr;
         _next_event = event;

      }
      _num_events++;
   }
}

void Custom_Sequence::add_event(Time_ms t, Event *event){

   Event *tmp = event->make_new_copy_of_event();
   tmp->_event_time = t;

   this->add_event(tmp);
}

/*****************************************************************************
 * Custom_Sequence::delete_event
 * Deletes the event with the given event_id from the event_list.
 *****************************************************************************/

int Custom_Sequence::delete_event(int event_id){

   Event *ptr = _find_event_id(event_id);
   if (ptr != NULL){
      _next_event->_next = ptr->_next;
      delete ptr;
      _num_events--;
      return 0;
   } else {
      return -1;
   }
}

/*****************************************************************************
 * Custom_Sequence::get_next_event_id
 * Returns the event_id of the next event in the queue.
 *****************************************************************************/

int Custom_Sequence::get_next_event_id(void){

   if (_next_event == NULL){
      return -1;
   } else {
      return _next_event->_event_id;
   }

}

/*****************************************************************************
 * Custom_Sequence::initialize_sequence
 * Initializes the pulse sequence and restores the spin model to equilibrium.
 *****************************************************************************/

Vector_3D& Custom_Sequence::initialize_sequence(Spin_Model& model){

   _next_event    = _event_list;
   model.restore_equilibrium();

   return (Vector_3D&)model;
}

/*****************************************************************************
 * Custom_Sequence::apply_to_time
 * Runs the pulse sequence up to the given time.
 *****************************************************************************/

Vector_3D& Custom_Sequence::apply_to_time(Spin_Model& model, Time_ms t){

#ifdef DEBUG
   assert(_next_event != NULL);
#endif

   while (t >= _next_event->_event_time){
      _next_event->apply(model);
      if (_next_event->_next == NULL){
         t -= _next_event->_event_time;
         _next_event = _event_list;
      } else {
         _next_event = _next_event->_next;
      }
   }
//   model.relax(t);
//   return (Vector_3D&)model;
   return model.get_time_sample(t);
   
}

/*****************************************************************************
 * Custom_Sequence::apply_next_event
 * Applies the next event in the pulse sequence.
 *****************************************************************************/

Vector_3D& Custom_Sequence::apply_next_event(Spin_Model& model){

#ifdef DEBUG
   assert(_next_event != NULL);
#endif

   _next_event->apply(model);
   _next_event = _next_event->_next;
   if (_next_event == NULL){
      _next_event = _event_list;
   }
   return model.get_net_magnetization();

}

/*****************************************************************************
 * Custom_Sequence::apply_to_end_of_repetition
 * Runs the pulse sequence up to the end of the repetition.
 *****************************************************************************/

Vector_3D& Custom_Sequence::apply_to_end_of_repetition(Spin_Model& model){

   // Apply to end of pulse sequence
   while(_next_event != NULL){
      _next_event->apply(model);
      _next_event = _next_event->_next;
   }
   // Once at end, return to beginning of pulse sequence
   _next_event = _event_list;

   return (Vector_3D&)model;
}

/*****************************************************************************
 * Custom_Sequence::apply_one_repetition
 * Runs the pulse sequence for one repetition.
 *****************************************************************************/

Vector_3D& Custom_Sequence::apply_one_repetition(Spin_Model& model){

   int n;
   for(n=0; n<_num_events; n++){
      _next_event->apply(model);
      _next_event = _next_event->_next;
      if (_next_event == NULL){
         _next_event = _event_list;
      }
   }

   return (Vector_3D&)model;
}

/*****************************************************************************
 * Custom_Sequence::apply_to_steady_state
 * Runs the pulse sequence up to steady state
 *****************************************************************************/

Vector_3D& Custom_Sequence::apply_to_steady_state(Spin_Model& model){

   Vector_3D steady_state;
   double    signal1, signal2;

   const double epsilon = 1E-4;

   steady_state = this->apply_one_repetition(model);
   signal1 = abs(steady_state);
   steady_state = this->apply_one_repetition(model);
   signal2 = abs(steady_state);

   while(fabs(signal1-signal2) > epsilon){
      signal1 = signal2;
      steady_state = this->apply_one_repetition(model);
      signal2 = abs(steady_state);
   }

   return (Vector_3D&)model;
}

/*****************************************************************************
 * Custom_Sequence::apply_trace
 * Trace the NMR signal response to a pulse sequence.
 *****************************************************************************/

void Custom_Sequence::apply_trace(Spin_Model& model, int trace_length, Time_ms trace_step,
                 Time_ms time[], Vector_3D m[]){

   int n;
   for(n=0; n<trace_length; n++){
      time[n] = n*trace_step;
      m[n]    = this->apply_to_time(model, time[n]);
   }
 
}

/*****************************************************************************
 * Custom_Sequence::dump_sequence_info
 * Dumps sequence information to a given file stream.
 *****************************************************************************/

void Custom_Sequence::display_sequence_info(ostream& stream){
   
   Event *ptr;
   char descriptor[80];

   this->display_info(stream);
   
   stream << "Programming Information:" << endl;
   stream << "------------------------" << endl;

   for (ptr=_event_list; ptr != NULL; ptr = ptr->_next){
      ptr->get_descriptor_string(descriptor);
      stream << descriptor << endl;
   }
}

   
/*****************************************************************************
 * Custom_Sequence protected member functions.
 *****************************************************************************/

/*****************************************************************************
 * Custom_Sequence::_delete_list
 * Deletes the list starting at position head.
 *****************************************************************************/

void Custom_Sequence::_delete_list(Event *head){

   Event *ptr;
   _next_event = head;

   while (_next_event != NULL){
      ptr = _next_event->_next;
      delete _next_event;
      _next_event = ptr;
   }

   head = (Event *)NULL;
   _next_event = _event_list;
}

/*****************************************************************************
 * Custom_Sequence::_find_event_id
 * Searches the event list for the event with the given event_id.
 * Returns a pointer to the required event and sets _current_event to
 * the previous event in the list.  If an event is not found returns NULL.
 *****************************************************************************/

Event *Custom_Sequence::_find_event_id(int event_id){

   _next_event = _event_list;

   Event *ptr = _next_event->_next;
   while(ptr != NULL && _next_event->_event_id != event_id){
      _next_event = ptr;
      ptr = ptr->_next;
   }
   return (ptr == NULL ? (Event *)NULL : _next_event);
}

/*****************************************************************************
 * Custom_Sequence::_find_event_after_time
 * Searches the event list for the last event with the given time.
 *
 * If an event with the given time is found, a pointer is returned
 * to the first event after the required time, and _current_event is
 * points to the last event with the given time.
 * If the given time is less than the time of the first event in the
 * list, a NULL pointer is returned and _current_event is set to NULL.
 * If the given time is greater than any time in the list, a NULL
 * pointer is returned and _current_event points to the last event
 * in the list (non-NULL).
 *
 *****************************************************************************/

Event *Custom_Sequence::_find_event_after_time(Time_ms t){

#ifdef DEBUG
   assert(_next_event != NULL);
#endif

   // Required time is less than the first event
   if (t < _event_list->_event_time){
      _next_event = (Event *)NULL;
      return (Event *)NULL;
   }

   // If the required time is less than the current event time
   // begin search at the beginning of the list
   if (t < _next_event->_event_time)
      _next_event = _event_list;

   Event *ptr = _next_event->_next;
   
   while(ptr != NULL && _next_event->_event_time <= t){
      _next_event = ptr;
      ptr = ptr->_next;
   }
   return (ptr == NULL ? (Event *)NULL : _next_event);
}
   


