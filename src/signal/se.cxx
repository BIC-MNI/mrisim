/*****************************************************************************
 * 
 * SE.CXX
 *
 * SE Class
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/se.cxx,v 1.1 2003-05-30 16:43:14 bert Exp $
// $Log: se.cxx,v $
// Revision 1.1  2003-05-30 16:43:14  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:33:42  rkwan
// Release 2.5
//
// Revision 2.1  1995/08/08  14:46:27  rkwan
// Updated for new spin model and pulse sequence.
//

#include "se.h"

/*****************************************************************************
 * SE Class
 * Pre-calculated Spin Echo steady state pulse sequence.
 *****************************************************************************/

/*****************************************************************************
 * SE constructors
 *****************************************************************************/

SE::SE() : Quick_Sequence() {
   _TR = 0;
   _TE = 0;
}

SE::SE(Time_ms TR, Time_ms TE) : Quick_Sequence() {
   _TR = TR;
   _TE = TE;
}

SE::SE(Time_ms TR, Time_ms TE, void (*f)(Vector_3D&,void *), void *data) :
   Quick_Sequence(f,data) {
   _TR = TR;
   _TE = TE;
}

/*****************************************************************************
 * SE destructor
 *****************************************************************************/

SE::~SE(){
}

/*****************************************************************************
 * SE::apply
 * Calculate the steady state magnetization of the spin echo sequence.
 *****************************************************************************/

void SE::apply(Quick_Model& model){

   Vector_3D signal;
   double NH;
   Time_ms T1, T2;

   NH = model.get_NH();
   T1 = model.get_T1();
   T2 = model.get_T2();

   signal[X_AXIS] = 0;
   signal[Y_AXIS] = NH*(1-2*exp(-(_TR-_TE/2)/T1)+exp(-_TR/T1))*exp(-_TE/T2);
   signal[Z_AXIS] = 0;

   model.set_net_magnetization(signal);
   (void)_sample->apply(model);
}

/*****************************************************************************
 * SE::display_info
 * Outputs information about the pulse sequence
 *****************************************************************************/

void SE::display_info(ostream& stream) const {

   stream << "Spin Echo:" << endl;
   stream << "TR         = " << _TR << " ms " << endl;
   stream << "TE         = " << _TE << " ms " << endl;

   Pulse_Sequence::display_info(stream);
}



