/*****************************************************************************
 * 
 * RF_PULSE.CXX
 *
 * RF_Pulse Class
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/rf_pulse.cxx,v 1.1 2003-05-30 16:43:13 bert Exp $
// $Log: rf_pulse.cxx,v $
// Revision 1.1  2003-05-30 16:43:13  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.3  1996/01/12  16:01:07  rkwan
// Added flip error.
//
// Revision 2.2  1995/12/11  14:27:12  rkwan
// Updated for fast_iso_model.
//
// Revision 2.1  1995/08/08  14:44:52  rkwan
// Updated for new spin model and pulse sequence.
//

#include <stdio.h>
#include "rf_pulse.h"

/*****************************************************************************
 * RF_Pulse Class
 *****************************************************************************/

/*****************************************************************************
 * RF_Pulse constructors
 *****************************************************************************/

RF_Pulse::RF_Pulse() : Event() {
   _axis  = X_AXIS;
   _angle = 0;
   _event_type = RF_PULSE;
}

RF_Pulse::RF_Pulse(Axis axis, Degrees angle) : Event() {
   _axis  = axis;
   _angle = angle;
   _event_type = RF_PULSE;
}

RF_Pulse::RF_Pulse(Axis axis, Degrees angle, Time_ms t) : Event(t) {
   _axis  = axis;
   _angle = angle;
   _event_type = RF_PULSE;
}

RF_Pulse::RF_Pulse(const RF_Pulse& rf) : Event((Event &)rf) {
   _axis = rf._axis;
   _angle = rf._angle;
   _event_type = RF_PULSE;
}

/*****************************************************************************
 * RF_Pulse::apply
 * Apply an RF pulse to a magnetization vector.
 *****************************************************************************/

Vector_3D& RF_Pulse::apply(Spin_Model& m){

   // Allow magnetization to relax to the event time.
   m.update(_event_time);        

   // Apply the RF pulse to rotate the magnetization.
   m.rotate(_event_time, m.get_flip_error()*_angle, _axis);

   return m.get_net_magnetization();
}

void RF_Pulse::get_descriptor_string(char s[]){

    char axis[7];
    switch(_axis){
       case X_AXIS:
          strcpy(axis,"X_AXIS");
          break;
       case Y_AXIS:
          strcpy(axis,"Y_AXIS");
          break;
       case Z_AXIS:
          strcpy(axis,"Z_AXIS");
          break; }

    sprintf((char *)s,"RF_Pulse %8.2lf %s %6.2lf",
            (double)_event_time,axis,(double)_angle);

} 

Event *RF_Pulse::make_new_copy_of_event(void){
   return (Event *)new RF_Pulse(*this);
}

