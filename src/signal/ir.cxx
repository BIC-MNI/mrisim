/*****************************************************************************
 * 
 * IR.CXX
 *
 * IR Class
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/ir.cxx,v 1.1 2003-05-30 16:43:13 bert Exp $
// $Log: ir.cxx,v $
// Revision 1.1  2003-05-30 16:43:13  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:32:49  rkwan
// Release 2.5
//
// Revision 2.1  1995/08/08  14:36:59  rkwan
// Updated for new spin model and pulse sequence.
//

#include "ir.h"

/*****************************************************************************
 * IR Class
 * Pre-calculated Inversion Recovery steady state pulse sequence.
 *****************************************************************************/

/*****************************************************************************
 * IR constructors
 *****************************************************************************/

IR::IR() : Quick_Sequence() {
   _TR = 0;
   _TE = 0;
   _TI = 0;
   _image_type = REAL_IMAGE;
}

IR::IR(Time_ms TR, Time_ms TE, Time_ms TI) : Quick_Sequence() {
   _TR = TR;
   _TE = TE;
   _TI = TI;
   _image_type = REAL_IMAGE;
}

IR::IR(Time_ms TR, Time_ms TE, Time_ms TI,
       void (*f)(Vector_3D&,void *), void *data) : Quick_Sequence(f,data) {
   _TR = TR;
   _TE = TE;
   _TI = TI;
   _image_type = REAL_IMAGE;
}

/*****************************************************************************
 * IR destructor
 *****************************************************************************/

IR::~IR(){
}

/*****************************************************************************
 * IR::apply
 * Calculate the steady state magnetization of the sequence.
 *****************************************************************************/

void IR::apply(Quick_Model& model){

   Vector_3D signal;
   double NH;
   Time_ms T1, T2;

   NH = model.get_NH();
   T1 = model.get_T1();
   T2 = model.get_T2();

   signal[X_AXIS] = 0;
   signal[Y_AXIS] = NH*(1-2*exp(-_TI/T1)+2*exp(-(_TR-_TE/2)/T1)
                         -exp(-_TR/T1))*exp(-_TE/T2);
   signal[Z_AXIS] = 0;

   model.set_net_magnetization(signal);
   (void)_sample->apply(model);
}

/*****************************************************************************
 * IR::display_info
 * Outputs information about the pulse sequence
 *****************************************************************************/

void IR::display_info(ostream& stream) const {

   stream << "Inversion Recovery:" << endl;
   stream << "TR         = " << _TR << " ms " << endl;
   stream << "TE         = " << _TE << " ms " << endl;
   stream << "TI         = " << _TI << " ms " << endl << endl;

   Pulse_Sequence::display_info(stream);
}



