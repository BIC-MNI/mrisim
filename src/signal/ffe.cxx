/*****************************************************************************
 * 
 * FFE.CXX
 *
 * FFE Class
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/ffe.cxx,v 1.1 2003-05-30 16:43:12 bert Exp $
// $Log: ffe.cxx,v $
// Revision 1.1  2003-05-30 16:43:12  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.2  1995/12/22  20:19:36  rkwan
// FFE sequences now default to 3D scan mode.
//
// Revision 2.1  1995/08/08  14:34:37  rkwan
// Updated for new spin model and pulse sequence.
//

#include "ffe.h"

/*****************************************************************************
 * FFE Class
 * Pre-calculated Fast Field Echo steady state pulse sequence.
 * This is a container class and cannot be instantiated.
 *****************************************************************************/

/*****************************************************************************
 * FFE constructors
 *****************************************************************************/

FFE::FFE() : Quick_Sequence() {
   _TR = 0;
   _TE = 0;
   _flip_angle = 0;
   _scan_mode = SCAN_MODE_3D;
}

FFE::FFE(Time_ms TR, Time_ms TE, Degrees flip_angle) : Quick_Sequence() {
   _TR = TR;
   _TE = TE;
   _flip_angle = flip_angle;
   _scan_mode = SCAN_MODE_3D;
}

FFE::FFE(Time_ms TR, Time_ms TE, Degrees flip_angle,
       void (*f)(Vector_3D&,void *), void *data) : Quick_Sequence(f,data) {
   _TR = TR;
   _TE = TE;
   _flip_angle = flip_angle;
   _scan_mode = SCAN_MODE_3D;
}

/*****************************************************************************
 * FFE destructor
 *****************************************************************************/

FFE::~FFE(){
}
