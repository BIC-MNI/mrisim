/*****************************************************************************
 * 
 * SPOILER.CXX
 *
 * Spoiler Class
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/spoiler.cxx,v 1.1 2003-05-30 16:43:14 bert Exp $
// $Log: spoiler.cxx,v $
// Revision 1.1  2003-05-30 16:43:14  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.2  1995/12/11  14:28:28  rkwan
// Updated for fast_iso_model.
//
// Revision 2.1  1995/08/08  14:54:42  rkwan
// Updated for new spin model and pulse sequence.
//

#include <stdio.h>
#include "spoiler.h"

/*****************************************************************************
 * Spoiler Class
 *****************************************************************************/

/*****************************************************************************
 * Spoiler::apply
 * Apply a spoiler (crusher) gradient a magnetization vector to
 * destroy all transverse magnetization.
 *****************************************************************************/

Vector_3D& Spoiler::apply(Spin_Model& m){

   // Allow magnetization to relax to the event time
   m.update(_event_time);

   // Destroy transverse magnetization
   m.zero_transverse_magnetization(_event_time);

   return m.get_net_magnetization();
}

void Spoiler::get_descriptor_string(char s[]){

   sprintf((char *)s,"Spoiler %8.2lf", (double)_event_time);

}

Event *Spoiler::make_new_copy_of_event(void){
   return (Event *)new Spoiler(*this);
}

