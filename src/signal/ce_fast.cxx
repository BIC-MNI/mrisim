/*****************************************************************************
 * 
 * CE_FAST.CXX
 *
 * CE_FAST Class
 *
 * R. Kwan
 * February 27, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/ce_fast.cxx,v 1.1 2003-05-30 16:43:12 bert Exp $
// $Log: ce_fast.cxx,v $
// Revision 1.1  2003-05-30 16:43:12  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:31:24  rkwan
// Release 2.5
//
// Revision 2.3  1996/01/12  16:00:57  rkwan
// Added flip error.
//
// Revision 2.2  1995/08/08  14:29:24  rkwan
// Updated for new spin model and pulse sequence.
//

#include "ce_fast.h"

/*****************************************************************************
 * CE_FAST Class
 * Pre-calculated CE_FAST steady state pulse sequence.
 *****************************************************************************/

/*****************************************************************************
 * CE_FAST::apply
 * Calculate the steady state magnetization of the sequence.
 *****************************************************************************/

void CE_FAST::apply(Quick_Model& model){

   Vector_3D signal;

   double  NH = model.get_NH();
   Time_ms T1 = model.get_T1();
   Time_ms T2 = model.get_T2();
   Time_ms T2s = model.get_T2s();

   double  E1 = exp(-_TR/T1);
   double  E2 = exp(-_TR/T2);
   double  cosa = cos(model.get_flip_error()*DEG_TO_RAD(_flip_angle));
   double  sina = sin(model.get_flip_error()*DEG_TO_RAD(_flip_angle));

   double  a = (1 - E1*E2*E2 + cosa*(E2*E2-E1))/(1-E1);
   double  b = (1 + cosa)*E2;

   signal[X_AXIS] = 0;
   signal[Y_AXIS] = (NH*sina/(1+cosa))*(1 + ((1+cosa)*E2*E2-a)/sqrt(a*a-b*b))*
                    exp((_TR-_TE)/T2s);
   signal[Z_AXIS] = 0;

   model.set_net_magnetization(signal);
   (void)_sample->apply(model);
}


/*****************************************************************************
 * CE_FAST::display_info
 * Outputs information about the pulse sequence
 *****************************************************************************/

void CE_FAST::display_info(ostream& stream) const {

   stream << "CE_FAST:" << endl;
   stream << "TR         = " << _TR << " ms " << endl;
   stream << "TE         = " << _TE << " ms " << endl;
   stream << "Flip angle = " << _flip_angle << " deg " << endl << endl;

   Pulse_Sequence::display_info(stream);
}


