/*****************************************************************************
 * 
 * SPOILED_FLASH.CXX
 *
 * Spoiled FLASH Class
 *
 * R. Kwan
 * February 27, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/spoiled_flash.cxx,v 1.1 2003-05-30 16:43:14 bert Exp $
// $Log: spoiled_flash.cxx,v $
// Revision 1.1  2003-05-30 16:43:14  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:33:57  rkwan
// Release 2.5
//
// Revision 2.2  1996/01/12  16:00:35  rkwan
// Added flip error.
//
// Revision 2.1  1995/08/08  14:52:08  rkwan
// Updated for new spin model and pulse sequence.
//

#include "spoiled_flash.h"

/*****************************************************************************
 * Spoiled FLASH Class
 * Pre-calculated Spoiled FLASH steady state pulse sequence.
 *****************************************************************************/

/*****************************************************************************
 * Spoiled_FLASH::apply
 * Calculate the steady state magnetization of the sequence.
 *****************************************************************************/

void Spoiled_FLASH::apply(Quick_Model& model){

   Vector_3D signal;

   double  NH = model.get_NH();
   Time_ms T1 = model.get_T1();
   Time_ms T2 = model.get_T2();
   Time_ms T2s = model.get_T2s();

   double  E1 = exp(-_TR/T1);
   double  flip_angle = model.get_flip_error()*DEG_TO_RAD(_flip_angle);

   signal[X_AXIS] = 0;
   signal[Y_AXIS] = NH*sin(flip_angle)*(1-E1)*exp(-_TE/T2s)/
                    (1 - E1*cos(flip_angle));
   signal[Z_AXIS] = 0;

   model.set_net_magnetization(signal);
   (void)_sample->apply(model);
}

/*****************************************************************************
 * Spoiled_FLASH::display_info
 * Outputs information about the pulse sequence
 *****************************************************************************/

void Spoiled_FLASH::display_info(ostream& stream) const {

   stream << "Spoiled_FLASH:" << endl;
   stream << "TR         = " << _TR << " ms " << endl;
   stream << "TE         = " << _TE << " ms " << endl;
   stream << "Flip angle = " << _flip_angle << " deg " << endl << endl;

   Pulse_Sequence::display_info(stream);
}




