/*****************************************************************************
 * 
 * SPIN_MODEL.CXX
 *
 * Container class for MRI simulation spin models.
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/spin_model.cxx,v 1.1 2003-05-30 16:43:14 bert Exp $
// $Log: spin_model.cxx,v $
// Revision 1.1  2003-05-30 16:43:14  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.4  1996/01/12  16:01:26  rkwan
// Added flip error.
//
// Revision 2.3  1995/12/11  15:07:30  rkwan
// Updated for fast_iso_model.
//
// Revision 2.2  1995/08/08  18:47:03  rkwan
// Added get_net_magnetization(float v[]) access function.
//
// Revision 2.1  1995/08/08  14:48:37  rkwan
// Updated for new spin model and pulse sequence.
//

#include "spin_model.h"

/*****************************************************************************
 * Spin_Model Class
 *****************************************************************************/

Spin_Model::Spin_Model() {
   _flip_error = 1.0;
}

Spin_Model::~Spin_Model() {}

Vector_3D& Spin_Model::get_net_magnetization(void) {
   return _net_magnetization; 
}

void Spin_Model::get_net_magnetization(float v[]){
   v[0] = _net_magnetization[X_AXIS];
   v[1] = _net_magnetization[Y_AXIS];
   v[2] = _net_magnetization[Z_AXIS];
}

Spin_Model::operator Vector_3D&() {
   return _net_magnetization;
}
