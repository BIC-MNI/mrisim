/*****************************************************************************
 * 
 * QUICK_MODEL.CXX
 *
 * Quick spin model for use with precalculated pulse sequences.
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/quick_model.cxx,v 1.1 2003-05-30 16:43:13 bert Exp $
// $Log: quick_model.cxx,v $
// Revision 1.1  2003-05-30 16:43:13  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.1  1995/08/08  14:41:43  rkwan
// Updated for new spin model and pulse sequence.
//

#include "quick_model.h"

/*****************************************************************************
 * Quick_Model Class
 *****************************************************************************/

/*****************************************************************************
 * Quick_Model Constructors
 *****************************************************************************/

Quick_Model::Quick_Model() : Spin_Model(), Tissue() {}

Quick_Model::Quick_Model(Time_ms T1, Time_ms T2, Time_ms T2s, float NH) :
   Spin_Model(), Tissue(T1,T2,T2s,NH) {}

Quick_Model::Quick_Model(const Tissue& tissue) :
   Spin_Model(), Tissue(tissue) {}

/*****************************************************************************
 * Quick_Model::set_net_magnetization
 *****************************************************************************/

void Quick_Model::set_net_magnetization(Vector_3D& v){
   _net_magnetization = v;
}



