/*****************************************************************************
 * 
 * REPEAT.CXX
 *
 * Repeat Class
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/repeat.cxx,v 1.1 2003-05-30 16:43:13 bert Exp $
// $Log: repeat.cxx,v $
// Revision 1.1  2003-05-30 16:43:13  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.2  1995/12/11  14:27:00  rkwan
// Updated for fast_iso_model.
//
// Revision 2.1  1995/08/08  14:44:13  rkwan
// Updated for new spin model and pulse sequence.
//

#include <stdio.h>
#include "repeat.h"

/*****************************************************************************
 * Repeat Class
 *****************************************************************************/

/*****************************************************************************
 * Repeat::apply
 * Repeat a pulse sequence.
 *****************************************************************************/

Vector_3D& Repeat::apply(Spin_Model& m){

   // Allow magnetization to relax to the event time
   m.update(_event_time);
   m.set_time(0.0);
   return m.get_net_magnetization();
}

void Repeat::get_descriptor_string(char s[]){

   sprintf((char *)s, "Repeat %8.2lf", (double)_event_time);

}

Event *Repeat::make_new_copy_of_event(void){
   return (Event *)new Repeat(*this);
}

