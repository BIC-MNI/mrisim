/*****************************************************************************
 * 
 * QUICKSEQ.CXX
 *
 * Quick_Sequence Class
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/quickseq.cxx,v 1.1 2003-05-30 16:43:13 bert Exp $
// $Log: quickseq.cxx,v $
// Revision 1.1  2003-05-30 16:43:13  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.1  1995/08/08  14:43:15  rkwan
// Updated for new spin model and pulse sequence.
//

#include "quickseq.h"

/*****************************************************************************
 * Quick_Sequence Class
 *****************************************************************************/

/*****************************************************************************
 * Quick_Sequence constructors
 *****************************************************************************/

Quick_Sequence::Quick_Sequence() : Pulse_Sequence() {
   _sample = new Sample;
}

Quick_Sequence::Quick_Sequence(void (*f)(Vector_3D&,void *), void *data) :
   Pulse_Sequence() {

   _sample = new Sample(0, f, data);

}

Quick_Sequence::Quick_Sequence(const Quick_Sequence& p) :
   Pulse_Sequence(p) {

   _sample = new Sample(*p._sample);

}

Quick_Sequence::~Quick_Sequence(){

   delete _sample;

}

