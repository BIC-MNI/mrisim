#ifndef __RF_PULSE_H
#define __RF_PULSE_H

/*****************************************************************************
 * 
 * RF_PULSE.H
 *
 * RF_Pulse Class
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/rf_pulse.h,v 1.1 2003-05-30 16:43:13 bert Exp $
 * $Log: rf_pulse.h,v $
 * Revision 1.1  2003-05-30 16:43:13  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.2  1995/12/11  14:27:07  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.1  1995/08/08  14:44:46  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "event.h"

/*****************************************************************************
 * RF_Pulse Class
 *****************************************************************************/

class RF_Pulse : public Event {
   public:
      RF_Pulse();                          // constructors
      RF_Pulse(Axis axis, Degrees angle);
      RF_Pulse(Axis axis, Degrees angle, Time_ms t);
      RF_Pulse(const RF_Pulse& rf);

      virtual Vector_3D& apply(Spin_Model& m);
      virtual void get_descriptor_string(char s[]);
      virtual Event *make_new_copy_of_event(void);

   protected:
      Axis    _axis;                       // apply RF pulse along this axis
      Degrees _angle;                      // rotation in degrees
};

#endif

