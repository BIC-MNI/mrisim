/*****************************************************************************
 * 
 * ISOCHROMAT_MODEL.CXX
 *
 * Isochromat spin system model class.
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/isochromat_model.cxx,v 1.1 2003-05-30 16:43:13 bert Exp $
// $Log: isochromat_model.cxx,v $
// Revision 1.1  2003-05-30 16:43:13  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:33:02  rkwan
// Release 2.5
//
// Revision 2.3  1996/01/17  17:58:22  rkwan
// Consistency fix between Fast_Isochromat_Model and Isochromat_Model.
//
// Revision 2.2  1995/12/11  14:24:02  rkwan
// Updated for fast_iso_model.
//
// Revision 2.1  1995/08/08  14:38:55  rkwan
// Updated for new spin model and pulse sequence.
//

#include "isochromat_model.h"

static const double MAGNITUDE_SCALING = 100000.0;

/*****************************************************************************
 * Isochromat_Model Class
 *****************************************************************************/

/*****************************************************************************
 * Isochromat_Model Constructors
 *****************************************************************************/

Isochromat_Model::Isochromat_Model(int n) : Spin_Model(), Tissue() {
   _num_of_isochromats = n;

   _m_equil            = new double[_num_of_isochromats];
   _off_resonance_freq = new Hertz[_num_of_isochromats];
   _m_0                = new Vector_3D[_num_of_isochromats];
   _m                  = new Vector_3D[_num_of_isochromats];

   _t_0     = (Time_ms)0.0;

   _net_magnetization = 0.0;

}

Isochromat_Model::Isochromat_Model(Time_ms T1, Time_ms T2, 
                                     Time_ms T2s, float NH, int n) :
   Spin_Model(), Tissue(T1,T2,T2s,NH) {

   // If number of isochromats is not specified use the default
   // exponential decay model
   if (n==0){
#ifdef DEBUG
       assert(_T2 != _T2s);
       assert(_T2s != 0.0);
#endif
       _bandwidth = (Hertz)(10*(_T2-_T2s)/(2*M_PI*_T2*_T2s));
       _num_of_isochromats = _next_power_of_two(
            (int)floor(40*_T2*_bandwidth));
   } else {
       _num_of_isochromats = n;
       _bandwidth = _num_of_isochromats/(40*_T2);
   }

   _m_equil            = new double[_num_of_isochromats];
   _off_resonance_freq = new Hertz[_num_of_isochromats];
   _m_0                = new Vector_3D[_num_of_isochromats];
   _m                  = new Vector_3D[_num_of_isochromats];

   _t_0     = (Time_ms)0.0;

   if (n==0){
      this->use_exponential_decay_model();
   }
   _net_magnetization = 0.0;

}

Isochromat_Model::Isochromat_Model(const Tissue& tissue, int n) :
   Spin_Model(), Tissue(tissue) {

   // If number of isochromats is not specified used the default
   // exponential decay model 
   if (n==0){
#ifdef DEBUG
       assert(_T2 != _T2s);
       assert(_T2s != 0.0);
#endif
       _bandwidth = (Hertz)(10*(_T2-_T2s)/(2*M_PI*_T2*_T2s));
       _num_of_isochromats = _next_power_of_two(
            (int)floor(40*_T2*_bandwidth));
   } else {
       _num_of_isochromats = n;
       _bandwidth = _num_of_isochromats/(40*_T2);
   }

#ifdef DEBUG
   cerr << "Isochromat_Model::_num_of_isochromats = " << _num_of_isochromats 
        << endl;
#endif

   _m_equil            = new double[_num_of_isochromats];
   _off_resonance_freq = new Hertz[_num_of_isochromats];
   _m_0                = new Vector_3D[_num_of_isochromats];
   _m                  = new Vector_3D[_num_of_isochromats];

   _t_0     = (Time_ms)0.0;

   if (n==0){
      this->use_exponential_decay_model();
   }
   _net_magnetization = 0.0;

}

/*****************************************************************************
 * Isochromat_Model Destructor
 *****************************************************************************/

Isochromat_Model::~Isochromat_Model(){
   delete[] _m_equil;
   delete[] _off_resonance_freq;
   delete[] _m_0;
   delete[] _m;
}

/*****************************************************************************
 * Isochromat_Model::restore_equilibrium
 *****************************************************************************/

void Isochromat_Model::restore_equilibrium(void){

   int n;

   _t_0 = (Time_ms)0.0;
   _net_magnetization = 0.0;

   for(n=0; n<_num_of_isochromats; n++){
      _m[n][X_AXIS] = 0.0;
      _m[n][Y_AXIS] = 0.0;
      _m[n][Z_AXIS] = _m_equil[n];
      _m_0[n] = _m[n];
      _net_magnetization += _m[n];
   }
   _net_magnetization *= (1.0/MAGNITUDE_SCALING);

}

/*****************************************************************************
 * Isochromat_Model::rotate
 *****************************************************************************/

void Isochromat_Model::rotate(Time_ms t, Degrees angle, Axis axis){
   int n;

   _t_0 = t;
   _net_magnetization = 0.0;
   for(n=0; n<_num_of_isochromats; n++){
      _m[n].rotate(axis, angle);
      _m_0[n] = _m[n];
      _net_magnetization += _m[n];
   }
   _net_magnetization *= (1.0/MAGNITUDE_SCALING);

}

void Isochromat_Model::rotate(Time_ms t, Degrees angle, Vector_3D& axis){
   int n;

   _t_0 = t;
   _net_magnetization = 0.0;
   for(n=0; n<_num_of_isochromats; n++){
      _m[n].rotate(axis, angle);
      _m_0[n] = _m[n];
      _net_magnetization += _m[n];
   }
   _net_magnetization *= (1.0/MAGNITUDE_SCALING);

}

/*****************************************************************************
 * Isochromat_Model::relax
 *****************************************************************************/

void Isochromat_Model::relax(Time_ms t){
   double E1, E2;

#ifdef DEBUG
   assert(t >= _t_0);
#endif

   E1 = exp(-(t-_t_0)/_T1);
   E2 = exp(-(t-_t_0)/_T2);

   int n;
   Vector_3D tmp;

   _net_magnetization = 0.0;

   for(n=0; n<_num_of_isochromats; n++){
      tmp = _m_0[n];
      tmp.rotate(Z_AXIS, 
                 RAD_TO_DEG(2*M_PI*(double)_off_resonance_freq[n]*(t-_t_0)));
      _m[n][X_AXIS] = E2*tmp[X_AXIS];
      _m[n][Y_AXIS] = E2*tmp[Y_AXIS];
      _m[n][Z_AXIS] = E1*tmp[Z_AXIS] + _m_equil[n]*(1-E1);
      _net_magnetization += _m[n];
   }

   _t = t;
   _net_magnetization *= (1.0/MAGNITUDE_SCALING);
}

void Isochromat_Model::update(Time_ms t){
   this->relax(t);
   int n;
   for(n=0; n<_num_of_isochromats; n++){
      _m_0[n] = _m[n];
   }
   _t_0 = t;
}
/*****************************************************************************
 * Isochromat_Model::use_linear_resonant_freq
 *****************************************************************************/

void Isochromat_Model::use_linear_resonant_freq(Hertz bandwidth){

   _bandwidth = bandwidth;
   Hertz Fs = (2*bandwidth)/_num_of_isochromats;
   unsigned int n;
   for(n=0; n<_num_of_isochromats; n++){
      _off_resonance_freq[n] = -bandwidth+n*Fs;
   }
   
}

/*****************************************************************************
 * Isochromat_Model::use_lorentzian_distribution
 *****************************************************************************/

void Isochromat_Model::use_lorentzian_distribution(double alpha){

#ifdef DEBUG
   assert(alpha != 0.0);
#endif

   unsigned int n;
   double Fs = 2*_bandwidth/_num_of_isochromats;

   for(n=0; n<_num_of_isochromats; n++){
      _m_equil[n] = 2*alpha*Fs*_NH*MAGNITUDE_SCALING/
                 (1.0+SQR(2*M_PI*alpha*(double)_off_resonance_freq[n]));
   }

}

/*****************************************************************************
 * Isochromat_Model::use_uniform_distribution
 *****************************************************************************/

void Isochromat_Model::use_uniform_distribution(void){
   int n;

   double Fs = 2*_bandwidth/_num_of_isochromats;
   for(n=0; n<_num_of_isochromats; n++){
      _m_equil[n] = Fs*MAGNITUDE_SCALING*_NH/_num_of_isochromats;
   }

}
   
/*****************************************************************************
 * Isochromat_Model::use_exponential_decay_model
 *****************************************************************************/

void Isochromat_Model::use_exponential_decay_model(void){
   double T21;
//   Hertz  offreson;

   T21 = _T2*_T2s/(_T2-_T2s);
//   offreson = (Hertz)((double)_num_of_isochromats/2.0)/(5*_T2);

   this->use_linear_resonant_freq(_bandwidth);
   this->use_lorentzian_distribution(T21);
}

/*****************************************************************************
 * Isochromat_Model::use_custom_resonant_freq
 *****************************************************************************/

void Isochromat_Model::use_custom_resonant_freq(double freq[]){

#ifdef DEBUG
   assert(freq != NULL);
#endif

   unsigned int n;
   for(n=0; n<_num_of_isochromats; n++){
      _off_resonance_freq[n] = (Hertz)freq[n];
   }

}

/*****************************************************************************
 * Isochromat_Model::use_custom_distribution
 *****************************************************************************/

void Isochromat_Model::use_custom_distribution(double dist[]){

#ifdef DEBUG
   assert(dist != NULL);
#endif

   unsigned int n;
   double norm_factor = 0.0;

   for(n=0; n<_num_of_isochromats; n++){
      _m_equil[n] = MAGNITUDE_SCALING*dist[n];
      norm_factor += _m_equil[n];
   }
   for(n=0; n<_num_of_isochromats; n++){
      _m_equil[n] = MAGNITUDE_SCALING*_NH*_m_equil[n]/norm_factor;
   }

}

/*****************************************************************************
 * Isochromat_Model::set_time
 *****************************************************************************/

void Isochromat_Model::set_time(Time_ms t){
   _t_0 = t;
}

/*****************************************************************************
 * Isochromat_Model::zero_transverse_magnetization
 *****************************************************************************/

void Isochromat_Model::zero_transverse_magnetization(Time_ms t){
   int n;

   for(n=0; n<_num_of_isochromats; n++){
      _m[n][X_AXIS] = 0.0;
      _m[n][Y_AXIS] = 0.0;
   }
   _t_0 = t;
}

/*****************************************************************************
 * Isochromat_Model::get_magnetization
 *****************************************************************************/

void Isochromat_Model::get_magnetization(int n, float v[]){

#ifdef DEBUG
   assert(n >= 0);
   assert(n <_num_of_isochromats);
#endif

   v[0] = _m[n][X_AXIS];
   v[1] = _m[n][Y_AXIS];
   v[2] = _m[n][Z_AXIS];

}

void Isochromat_Model::display_model_info(ostream& stream){

   stream << "Isochromat Spin Model Info:" << endl;
   stream << "Number of isochromats: " << _num_of_isochromats << endl;
   stream << "T1:   " << _T1 << endl;
   stream << "T2:   " << _T2 << endl;
   stream << "T2*:  " << _T2s << endl;
   stream << "NH:   " << _NH << endl;

}

void Isochromat_Model::display_model_state(ostream& stream){
   stream << "Isochromat Spin Model State:" << endl;
   stream << "_num_of_isochromats: " << _num_of_isochromats << endl;
   stream << "_t_0: " << _t_0 << " _t: " << _t << endl;
   stream << "_net_magnetization: ";
   _net_magnetization.print(stream);
}
