/*****************************************************************************
 * 
 * VECTOR_MODEL.CXX
 *
 * Single magnetization vector model.
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/vector_model.cxx,v 1.1 2003-05-30 16:43:14 bert Exp $
// $Log: vector_model.cxx,v $
// Revision 1.1  2003-05-30 16:43:14  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 1.1  1995/12/11  14:32:20  rkwan
// Initial revision
//

#include "vector_model.h"

/*****************************************************************************
 * Vector_Model Class
 *****************************************************************************/

/*****************************************************************************
 * Vector_Model Constructors
 *****************************************************************************/

Vector_Model::Vector_Model() : Spin_Model(), Tissue() {
   _m_equil = (double)_NH;
   _t_0     = (Time_ms)0.0;
   this->restore_equilibrium();
}

Vector_Model::Vector_Model(Time_ms T1, Time_ms T2, 
                                     Time_ms T2s, float NH) :
   Spin_Model(), Tissue(T1,T2,T2s,NH) {
   
   _m_equil = (double)_NH;
   _t_0     = (Time_ms)0.0;
   this->restore_equilibrium();
}

Vector_Model::Vector_Model(const Tissue& tissue) :
   Spin_Model(), Tissue(tissue) {
   
   _m_equil = (double)_NH;
   _t_0     = (Time_ms)0.0;
   this->restore_equilibrium();
}

/*****************************************************************************
 * Vector_Model::restore_equilibrium
 *****************************************************************************/

void Vector_Model::restore_equilibrium(void){
   _net_magnetization[X_AXIS] = 0.0;
   _net_magnetization[Y_AXIS] = 0.0;
   _net_magnetization[Z_AXIS] = _m_equil;
   _m_0 = _net_magnetization;
   _t_0       = (Time_ms)0.0;
}

void Vector_Model::rotate(Time_ms t, Degrees angle, Axis axis){
   _net_magnetization.rotate(axis, angle);
   _t_0 = t;
   _m_0 = _net_magnetization;
}

void Vector_Model::rotate(Time_ms t, Degrees angle, Vector_3D& axis){
   _net_magnetization.rotate(axis, angle);
   _t_0 = t;
   _m_0 = _net_magnetization;
}

void Vector_Model::relax(Time_ms t){
   double E1, E2;

#ifdef DEBUG
   assert(t >= _t_0);
#endif

   E1 = exp(-(t-_t_0)/_T1);
   E2 = exp(-(t-_t_0)/_T2);

   _net_magnetization[X_AXIS] = E2*_m_0[X_AXIS];
   _net_magnetization[Y_AXIS] = E2*_m_0[Y_AXIS];
   _net_magnetization[Z_AXIS] = E1*_m_0[Z_AXIS] + _m_equil*(1-E1);

   _t = t;
}

void Vector_Model::update(Time_ms t){
   this->relax(t);
   _m_0 = _net_magnetization;
   _t_0 = t;
}

void Vector_Model::set_time(Time_ms t){
   _t_0 = t;
}

void Vector_Model::zero_transverse_magnetization(Time_ms t){
   _net_magnetization[X_AXIS] = 0.0;
   _net_magnetization[Y_AXIS] = 0.0;
   _t_0 = t;
}

void Vector_Model::display_model_info(ostream& stream){

   stream << "Single Spin Model Info:" << endl;
   stream << "Name: " << _tissue_name << endl;
   stream << "T1:   " << _T1 << endl;
   stream << "T2:   " << _T2 << endl;
   stream << "T2*:  " << _T2s << endl;
   stream << "NH:   " << _NH << endl;

}

void Vector_Model::display_model_state(ostream& stream){
   stream << "Single Spin Model State:" << endl;
   stream << "_net_magnetization: ";
   _net_magnetization.print(stream);
   stream << "_m_equil: " << _m_equil << endl;
   stream << "_m_0: ";
   _m_0.print(stream);
   stream << "_t_0: " << _t_0 << " _t: " << _t << endl;
}
