/*****************************************************************************
 * 
 * FAST_ISO_MODEL.CXX
 *
 * Fast Isochromat spin system model uses FFT.
 *
 * R. Kwan
 * November 22, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/fast_iso_model.cxx,v 1.1 2003-05-30 16:43:12 bert Exp $
// $Log: fast_iso_model.cxx,v $
// Revision 1.1  2003-05-30 16:43:12  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:31:57  rkwan
// Release 2.5
//
// Revision 1.4  1996/01/17  17:57:48  rkwan
// Consistency fix between Fast_Isochromat_Model and Isochromat_Model.
//
// Revision 1.3  1995/12/22  20:19:23  rkwan
// Update for percent coil and RF map features.
//
// Revision 1.2  1995/12/11  20:58:04  rkwan
// Fix FFT/chirp relaxation.
//
// Revision 1.1  1995/12/11  14:23:38  rkwan
// Initial revision
//

#include "fast_iso_model.h"
#include <math.h>

static const double MAGNITUDE_SCALING = 100000.0;

/*****************************************************************************
 * Fast_Isochromat_Model Class
 *****************************************************************************/

/*****************************************************************************
 * Fast_Isochromat_Model Constructors
 *****************************************************************************/

Fast_Isochromat_Model::Fast_Isochromat_Model(int n) : Spin_Model(), Tissue() {
   _num_of_isochromats = n;

   _m_equil            = new double[_num_of_isochromats];
   _off_resonance_freq = new Hertz[_num_of_isochromats];

   _m_0_z              = new double[_num_of_isochromats];
   _m_0_xy             = new double[2*_num_of_isochromats];
   _m_z                = new double[_num_of_isochromats];
   _m_xy               = new double[2*_num_of_isochromats];

   _t_0     = (Time_ms)0.0;

}

Fast_Isochromat_Model::Fast_Isochromat_Model(Time_ms T1, Time_ms T2, 
                                     Time_ms T2s, float NH, int n) :
   Spin_Model(), Tissue(T1,T2,T2s,NH) {

   // If number of isochromats is not specified use the default
   // exponential decay model 
   if (n == 0){
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
   _m_0_z              = new double[_num_of_isochromats];
   _m_0_xy             = new double[2*_num_of_isochromats];
   _m_z                = new double[_num_of_isochromats];
   _m_xy               = new double[2*_num_of_isochromats];

   _t_0     = (Time_ms)0.0;

   this->use_exponential_decay_model();
    
}

Fast_Isochromat_Model::Fast_Isochromat_Model(const Tissue& tissue, int n) :
   Spin_Model(), Tissue(tissue) {
   
   // If number of isochromats is not specified used the default
   // exponential decay model 
   if (n == 0){
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
   _m_0_z              = new double[_num_of_isochromats];
   _m_0_xy             = new double[2*_num_of_isochromats];
   _m_z                = new double[_num_of_isochromats];
   _m_xy               = new double[2*_num_of_isochromats];

   _t_0     = (Time_ms)0.0;

   this->use_exponential_decay_model();

}

/*****************************************************************************
 * Fast_Isochromat_Model Destructor
 *****************************************************************************/

Fast_Isochromat_Model::~Fast_Isochromat_Model(){
   delete[] _m_equil;
   delete[] _off_resonance_freq;
   delete[] _m_0_z;
   delete[] _m_0_xy;
   delete[] _m_z;
   delete[] _m_xy;
}

/*****************************************************************************
 * Fast_Isochromat_Model::restore_equilibrium
 *****************************************************************************/

void Fast_Isochromat_Model::restore_equilibrium(void){

   // Clear transverse magnetization and restore z to equilibrium
   memset(_m_xy, 0, 2*_num_of_isochromats*sizeof(double));
   memcpy(_m_z, _m_equil, _num_of_isochromats*sizeof(double));

   // Restore initial magnetization
   _update_initial_mag((Time_ms)0.0);

   // Compute net magnetization
   _compute_net_mag(_net_magnetization, _m_z, _m_xy);

}

/*****************************************************************************
 * Fast_Isochromat_Model::rotate
 *****************************************************************************/

void Fast_Isochromat_Model::rotate(Time_ms t, Degrees angle, Axis axis){
   double m_x, m_y, m_z;
   double rotation[3][3];
   unsigned int n;

   // Compute rotation matrix
   _compute_rotation_matrix(rotation, axis, angle);
   
   // Apply rotation matrix to each isochromat   
   for(n=0; n<_num_of_isochromats; n++){
      m_y = _m_xy[2*n];
      m_x = _m_xy[2*n+1];
      m_z = _m_z[n];
      _m_xy[2*n+1] = rotation[0][0]*m_x + rotation[0][1]*m_y + 
                     rotation[0][2]*m_z; 
      _m_xy[2*n]   = rotation[1][0]*m_x + rotation[1][1]*m_y + 
                     rotation[1][2]*m_z; 
      _m_z[n]      = rotation[2][0]*m_x + rotation[2][1]*m_y + 
                     rotation[2][2]*m_z; 
   }

   _update_initial_mag(t);
   _compute_net_mag(_net_magnetization, _m_z, _m_xy);

}

void Fast_Isochromat_Model::rotate(Time_ms t, Degrees angle, Vector_3D& axis){
   double m_x, m_y, m_z;
   double rotation[3][3];
   unsigned int n;

   // Compute rotation matrix
   _compute_general_rotation_matrix(rotation, axis, angle);

   // Apply rotation matrix to each isochromat   
   for(n=0; n<_num_of_isochromats; n++){
      m_y = _m_xy[2*n];
      m_x = _m_xy[2*n+1];
      m_z = _m_z[n];
      _m_xy[2*n+1] = rotation[0][0]*m_x + rotation[0][1]*m_y + 
                     rotation[0][2]*m_z; 
      _m_xy[2*n]   = rotation[1][0]*m_x + rotation[1][1]*m_y + 
                     rotation[1][2]*m_z; 
      _m_z[n]      = rotation[2][0]*m_x + rotation[2][1]*m_y + 
                     rotation[2][2]*m_z; 
   }

   _update_initial_mag(t);
   _compute_net_mag(_net_magnetization, _m_z, _m_xy);

}

/*****************************************************************************
 * Fast_Isochromat_Model::relax
 *****************************************************************************/

void Fast_Isochromat_Model::relax(Time_ms t){

#ifdef DEBUG
   assert(t >= _t_0);
#endif

   double E1 = exp(-(t-_t_0)/_T1);
   double E2 = exp(-(t-_t_0)/_T2);

   double cosa, sina;
   unsigned int n;

   // Compute transverse magnetization
   // Faster trig recursion algorithm

   // Compute angle increment
   cosa = cos(2*M_PI*(2*_bandwidth/(double)_num_of_isochromats)*(t-_t_0));
   sina = sin(2*M_PI*(2*_bandwidth/(double)_num_of_isochromats)*(t-_t_0));

   // Compute initial cos/sin to begin recursion
   double c1 = cos(2*M_PI*(double)_off_resonance_freq[0]*(t-_t_0));
   double s1 = sin(2*M_PI*(double)_off_resonance_freq[0]*(t-_t_0));
   double c2, s2;

   for(n=0; n<_num_of_isochromats; n+=2){
      // Compute dephasing rotation and relaxation
      _m_xy[n+1] = E2*(_m_0_xy[n+1]*c1 - _m_0_xy[n]*s1);
      _m_xy[n]   = E2*(_m_0_xy[n+1]*s1 + _m_0_xy[n]*c1);

      // Update rotation using recursion relation
      c2 = cosa*c1 - sina*s1;
      s2 = sina*c1 + cosa*s1;
      c1 = c2;  
      s1 = s2;
   }

   // Account for phase inversion due to fftshift half way through
   // data vector.  cos(-x)=cos(x), sin(-x) = -sin(x)
 
   s1 = -s2;
   for(n=_num_of_isochromats; n<2*_num_of_isochromats; n+=2){
      // Compute dephasing rotation and relaxation
      _m_xy[n+1] = E2*(_m_0_xy[n+1]*c1 - _m_0_xy[n]*s1);
      _m_xy[n]   = E2*(_m_0_xy[n+1]*s1 + _m_0_xy[n]*c1);

      // Update rotation using recursion relation
      c2 = cosa*c1 - sina*s1;
      s2 = sina*c1 + cosa*s1;
      c1 = c2;  
      s1 = s2;
   }

   // Compute longitudinal magnetization
   for(n=0; n<_num_of_isochromats; n++){
      _m_z[n]      = E1*_m_0_z[n] + _m_equil[n]*(1-E1);
   }
 
   _compute_net_mag(_net_magnetization, _m_z, _m_xy);
   _t = t;
}

void Fast_Isochromat_Model::relax_slow(Time_ms t){

#ifdef DEBUG
   assert(t >= _t_0);
#endif

   double E1 = exp(-(t-_t_0)/_T1);
   double E2 = exp(-(t-_t_0)/_T2);

   double cosa, sina;
   unsigned int n;

   // Compute transverse magnetization
   for(n=0; n<_num_of_isochromats; n++){
      // Compute dephasing rotation
      cosa = cos(2*M_PI*(double)_off_resonance_freq[n]*(t-_t_0));
      sina = sin(2*M_PI*(double)_off_resonance_freq[n]*(t-_t_0));

      // Compute rotation and relaxation
      _m_xy[2*n+1] = E2*(_m_0_xy[2*n+1]*cosa - _m_0_xy[2*n]*sina);
      _m_xy[2*n]   = E2*(_m_0_xy[2*n+1]*sina + _m_0_xy[2*n]*cosa);
   }

   // Compute longitudinal magnetization
   for(n=0; n<_num_of_isochromats; n++){
      _m_z[n]      = E1*_m_0_z[n] + _m_equil[n]*(1-E1);
   }
 
   _compute_net_mag(_net_magnetization, _m_z, _m_xy);
   _t = t;
}

/*****************************************************************************
 * Fast_Isochromat_Model::update
 *****************************************************************************/

void Fast_Isochromat_Model::update(Time_ms t){
   this->relax(t);
   this->_update_initial_mag(t);
}

/*****************************************************************************
 * Fast_Isochromat_Model::use_linear_resonant_freq
 *****************************************************************************/

void Fast_Isochromat_Model::use_linear_resonant_freq(Hertz bandwidth){

   _bandwidth = bandwidth;
   Hertz Fs = (2*bandwidth)/_num_of_isochromats;
   unsigned int n;
   for(n=0; n<_num_of_isochromats; n++){
      _off_resonance_freq[n] = -bandwidth+n*Fs;
   }
   _fftshift_1d(_off_resonance_freq, _num_of_isochromats, sizeof(double));

}

/*****************************************************************************
 * Fast_Isochromat_Model::use_lorentzian_distribution
 *****************************************************************************/

void Fast_Isochromat_Model::use_lorentzian_distribution(double alpha){

#ifdef DEBUG
   assert(alpha != 0.0);
#endif

   unsigned int n;
   double Fs = 2*_bandwidth/_num_of_isochromats;

   for(n=0; n<_num_of_isochromats; n++){
      _m_equil[n] = 2*alpha*_NH*Fs*MAGNITUDE_SCALING/
                  (1.0+SQR(2*M_PI*alpha*(double)_off_resonance_freq[n]));
   }

}

/*****************************************************************************
 * Fast_Isochromat_Model::use_uniform_distribution
 *****************************************************************************/

void Fast_Isochromat_Model::use_uniform_distribution(void){
   unsigned int n;

   double Fs = 2*_bandwidth/_num_of_isochromats;
   for(n=0; n<_num_of_isochromats; n++){
      _m_equil[n] = Fs*MAGNITUDE_SCALING*_NH/_num_of_isochromats;
   }

}
   
/*****************************************************************************
 * Fast_Isochromat_Model::use_exponential_decay_model
 *****************************************************************************/

void Fast_Isochromat_Model::use_exponential_decay_model(void){
   double T21;
//   Hertz  offreson;

   T21 = _T2*_T2s/(_T2-_T2s);
//   offreson = (Hertz)((double)_num_of_isochromats/2.0)/(5*_T2);

   this->use_linear_resonant_freq(_bandwidth);
   this->use_lorentzian_distribution(T21);
}

/*****************************************************************************
 * Fast_Isochromat_Model::use_custom_resonant_freq
 *****************************************************************************/

void Fast_Isochromat_Model::use_custom_resonant_freq(double freq[]){

#ifdef DEBUG
   assert(freq != NULL);
#endif

   unsigned int n;
   for(n=0; n<_num_of_isochromats; n++){
      _off_resonance_freq[n] = (Hertz)freq[n];
   }
   _fftshift_1d(_off_resonance_freq, _num_of_isochromats, sizeof(double));

}

/*****************************************************************************
 * Fast_Isochromat_Model::use_custom_distribution
 *****************************************************************************/

void Fast_Isochromat_Model::use_custom_distribution(double dist[]){

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
 * Fast_Isochromat_Model::set_time
 *****************************************************************************/

void Fast_Isochromat_Model::set_time(Time_ms t){
   _t_0 = t;
}

/*****************************************************************************
 * Fast_Isochromat_Model::zero_transverse_magnetization
 *****************************************************************************/

void Fast_Isochromat_Model::zero_transverse_magnetization(Time_ms t){
   memset(_m_xy, 0, 2*_num_of_isochromats*sizeof(double));
   _t_0 = t;
}

/*****************************************************************************
 * Fast_Isochromat_Model::spoil
 *****************************************************************************/

void Fast_Isochromat_Model::spoil(double G, Time_ms t){
   unsigned int n;

   // Compute angle increment
   double cosa = cos(2*M_PI*2*_bandwidth*t*G/_num_of_isochromats);
   double sina = sin(2*M_PI*2*_bandwidth*t*G/_num_of_isochromats);

   // Compute initial cos/sin to begin recursion
   double c1 = cos(2*M_PI*_off_resonance_freq[0]*G*t);
   double s1 = sin(2*M_PI*_off_resonance_freq[0]*G*t);
   double c2, s2;

   for(n=0; n<_num_of_isochromats; n+=2){
      // Compute dephasing rotation 
      _m_xy[n+1] = (_m_0_xy[n+1]*c1 - _m_0_xy[n]*s1);
      _m_xy[n]   = (_m_0_xy[n+1]*s1 + _m_0_xy[n]*c1);

      // Update rotation using recursion relation
      c2 = cosa*c1 - sina*s1;
      s2 = sina*c1 + cosa*s1;
      c1 = c2;  
      s1 = s2;
   }

   // Account for phase inversion due to fftshift half way through
   // data vector.  cos(-x)=cos(x), sin(-x) = -sin(x)
    
   s1 = -s2;
   for(n=_num_of_isochromats; n<2*_num_of_isochromats; n+=2){
      // Compute dephasing rotation
      _m_xy[n+1] = (_m_0_xy[n+1]*c1 - _m_0_xy[n]*s1);
      _m_xy[n]   = (_m_0_xy[n+1]*s1 + _m_0_xy[n]*c1);
   
      // Update rotation using recursion relation
      c2 = cosa*c1 - sina*s1;
      s2 = sina*c1 + cosa*s1;
      c1 = c2;  
      s1 = s2;
   }

/* OBSOLETE CODE */
/*
   for(n=0; n<_num_of_isochromats; n++){
      cosa = cos((double)_off_resonance_freq[n]*G*t);
      sina = sin((double)_off_resonance_freq[n]*G*t);
      m_y = _m_xy[2*n];
      m_x = _m_xy[2*n+1]; 
      _m_xy[2*n]   = m_x*sina + m_y*cosa;
      _m_xy[2*n+1] = m_x*cosa - m_y*sina;
   }
*/
   _compute_net_mag(_net_magnetization, _m_z, _m_xy);
   _update_initial_mag(t);
}

/*****************************************************************************
 * Fast_Isochromat_Model::get_magnetization
 *****************************************************************************/

void Fast_Isochromat_Model::get_magnetization(int n, float v[]){

#ifdef DEBUG
   assert(n >= 0);
   assert(n <_num_of_isochromats);
#endif

   v[0] = _m_xy[2*n+1];
   v[1] = _m_xy[2*n];
   v[2] = _m_z[n];

}

/*****************************************************************************
 * Fast_Isochromat_Model::get_time_samples
 *****************************************************************************/

void Fast_Isochromat_Model::get_time_samples(double xy_samples[],
                                             double z_samples[]){

   // Time scale step is inverse of frequency bandwidth
   // Assume last entry is the highest resonance frequency
   double time_step = (double)this->get_time_step();

   // Copy initial magnetization
   memcpy(xy_samples, _m_0_xy, 2*_num_of_isochromats*sizeof(double));

   // Compute time samples using FFT
   four1(xy_samples-1, _num_of_isochromats, -1);

   // Compute z-magnetization
   unsigned int n;
   double m_equil = 0.0, m_0_z = 0.0;
   for(n=0; n<_num_of_isochromats; n++){
      m_0_z   += _m_0_z[n];
      m_equil += _m_equil[n];
   }

   // Compensate for relaxation

   double E1, E2;
   for(n=0; n<_num_of_isochromats; n++){
      E2 = exp(-(double)n*time_step/_T2)/MAGNITUDE_SCALING;
      xy_samples[2*n]   *= E2;
      xy_samples[2*n+1] *= E2;
   }

   for(n=0; n<_num_of_isochromats; n++){
      E1 = exp(-(double)n*time_step/_T1);
      z_samples[n]  = (E1*m_0_z + m_equil*(1-E1))/MAGNITUDE_SCALING;
   }

}

void Fast_Isochromat_Model::get_time_samples(double xy_samples[],
                                             double z_samples[],
                                             int num_of_samples,
                                             Time_ms t_start,
                                             Time_ms t_step){

   // Calculate Chirp parameters
   double freq_step = 2*_bandwidth/_num_of_isochromats;
   double start = 2*M_PI*t_start*freq_step;
   double increment = 2*M_PI*t_step*freq_step;

   // Compute time samples using Chirp-DFT algorithm
   chirpdft(_m_0_xy, _num_of_isochromats, xy_samples, num_of_samples,
            start, increment);

   // Compute z-magnetization
   unsigned int n;
   double m_equil = 0.0, m_0_z = 0.0;
   for(n=0; n<_num_of_isochromats; n++){
      m_0_z   += _m_0_z[n];
      m_equil += _m_equil[n];
   }

   // Compensate for relaxation
   double E1, E2;
   double scale = 1.0/MAGNITUDE_SCALING;

   for(n=0; n<num_of_samples; n++){
      E2 = exp(-(t_start+(double)n*t_step)/_T2)*scale;
      xy_samples[2*n]   *= E2;
      xy_samples[2*n+1] *= E2;
   }

   for(n=0; n<num_of_samples; n++){
      E1 = exp(-(t_start+(double)n*t_step)/_T1);
      z_samples[n]  = (E1*m_0_z + m_equil*(1-E1))*scale;
   }
}

/*****************************************************************************
 * Fast_Isochromat_Model::get_time_sample
 *****************************************************************************/

void Fast_Isochromat_Model::get_time_sample(Vector_3D& sample, 
                                            int sample_num){

   // Time scale step is inverse of frequency bandwidth
   // Assume last entry is the highest resonance frequency
   Time_ms time_step = fabs(this->get_time_step());

   // Compute time sample using Goertzel DFT algorithm
   double tmp[2];
   int    points[1] = {sample_num};
   dftgl(_m_0_xy, _num_of_isochromats, tmp, 1, points);

   // Compute z-magnetization
   unsigned int n;
   double m_equil = 0.0, m_0_z = 0.0;
   for(n=0; n<_num_of_isochromats; n++){
      m_0_z   += _m_0_z[n];
      m_equil += _m_equil[n];
   }

   // Compensate for relaxation
   double E1 = exp(-sample_num*time_step/_T1);
   double E2 = exp(-sample_num*time_step/_T2);
   sample[X_AXIS] = tmp[1]*E2/MAGNITUDE_SCALING;
   sample[Y_AXIS] = tmp[0]*E2/MAGNITUDE_SCALING;
   sample[Z_AXIS] = (E1*m_0_z + m_equil*(1-E1))/MAGNITUDE_SCALING;

}

Vector_3D& Fast_Isochromat_Model::get_time_sample(Time_ms t){

   // Time scale step is inverse of frequency bandwidth
   // Assume last entry is the highest resonance frequency
   Time_ms time_step = fabs(this->get_time_step());
   int k = (int)floor((t-_t_0)/time_step);

   // Compute time sample using Goertzel DFT algorithm
   double tmp[4];
   int    points[2] = {k, k+1};

   dftgl(_m_0_xy, _num_of_isochromats, tmp, 2, points);

   // Compute z-magnetization
   unsigned int n;
   double m_equil = 0.0, m_0_z = 0.0;
   for(n=0; n<_num_of_isochromats; n++){
      m_0_z   += _m_0_z[n];
      m_equil += _m_equil[n];
   }

   // Compensate for relaxation
   double E1 = exp(-(t-_t_0)/_T1);
   double E2 = exp(-(t-_t_0)/_T2);
   double scale = (t-_t_0-k*time_step)/time_step;
   _net_magnetization[X_AXIS] = E2*(scale*(tmp[3]-tmp[1])+tmp[1])/
                                MAGNITUDE_SCALING;
   _net_magnetization[Y_AXIS] = E2*(scale*(tmp[2]-tmp[0])+tmp[0])/
                                MAGNITUDE_SCALING;
   _net_magnetization[Z_AXIS] = (E1*m_0_z + m_equil*(1-E1))/
                                MAGNITUDE_SCALING;

   return _net_magnetization;

}

/*****************************************************************************
 * Fast_Isochromat_Model::display_model_info
 *****************************************************************************/

void Fast_Isochromat_Model::display_model_info(ostream& stream){

   stream << "Fast Isochromat Spin Model:" << endl;
   stream << "---------------------------" << endl;
   stream << "Number of Isochromats: " << _num_of_isochromats << endl;
   stream << "Bandwidth:             " << _bandwidth << endl;
   stream << "Time Step:             " << 1/_bandwidth << endl;
   stream << "T1:                    " << _T1 << endl;
   stream << "T2:                    " << _T2 << endl;
   stream << "T2*:                   " << _T2s << endl; 
   stream << "NH:                    " << _NH << endl;
   stream << endl;
}

/*****************************************************************************
 * Fast_Isochromat_Model::display_model_state
 *****************************************************************************/

void Fast_Isochromat_Model::display_model_state(ostream& stream){
   stream << "Fast Isochromat Spin Model State:" << endl;
   stream << "_num_of_isochromats: " << _num_of_isochromats << endl;
   stream << "_t_0: " << _t_0 << " _t: " << _t << endl;
   stream << "_net_magnetization: ";
   _net_magnetization.print(stream);
}

/*****************************************************************************
 * Fast_Isochromat_Model private member functions
 *****************************************************************************/

/*****************************************************************************
 * Fast_Isochromat_Model::_compute_net_mag
 *****************************************************************************/

void Fast_Isochromat_Model::_compute_net_mag(Vector_3D& net, 
                                             double z[], double xy[]){
   double m_x = 0.0;
   double m_y = 0.0;
   double m_z = 0.0;
   unsigned int n;

   for(n=0; n<_num_of_isochromats; n++){
      m_z += z[n];
   }
   for(n=0; n<2*_num_of_isochromats; n+=2){
      m_y += xy[n];
      m_x += xy[n+1];
   }
   net[X_AXIS] = m_x/MAGNITUDE_SCALING;
   net[Y_AXIS] = m_y/MAGNITUDE_SCALING;
   net[Z_AXIS] = m_z/MAGNITUDE_SCALING;
}

/*****************************************************************************
 * Fast_Isochromat_Model::_update_initial_mag
 *****************************************************************************/

void Fast_Isochromat_Model::_update_initial_mag(Time_ms t){

   memcpy(_m_0_z, _m_z, _num_of_isochromats*sizeof(double));
   memcpy(_m_0_xy, _m_xy, 2*_num_of_isochromats*sizeof(double));
   _t_0 = t;

}

/*****************************************************************************
 * Fast_Isochromat_Model::_compute_rotation_matrix
 *****************************************************************************/

void Fast_Isochromat_Model::_compute_rotation_matrix(double rotation[][3],
     Axis axis, Degrees angle){

   // Clear rotation matrix
   memset(rotation, 0, 9*sizeof(double));

   // Precompute sin/cos
   double sina, cosa;
   if (angle == 90){
      sina = 1; cosa = 0;
   } else if (angle == 180){
      sina = 0; cosa = -1;
   } else {
      sina = sin(DEG_TO_RAD(angle));
      cosa = cos(DEG_TO_RAD(angle));
   }

   // Compute rotation matrices
   switch(axis){
      case X_AXIS:
         rotation[0][0] = 1.0;
         rotation[1][1] = cosa;
         rotation[1][2] = sina;
         rotation[2][1] = -sina;
         rotation[2][2] = cosa;
         break;
      case Y_AXIS: 
         rotation[0][0] = cosa;
         rotation[0][2] = -sina;
         rotation[1][1] = 1.0;
         rotation[2][0] = sina;
         rotation[2][2] = cosa;
         break;
      case Z_AXIS:
         rotation[0][0] = cosa; 
         rotation[0][1] = sina;
         rotation[1][0] = -sina;
         rotation[1][1] = cosa;
         rotation[2][2] = 1.0;
         break;
   }
}

void Fast_Isochromat_Model::_compute_general_rotation_matrix(
     double rotation[][3],
     Vector_3D& axis, Degrees angle){

   // Make axis a unit vector
   axis[X_AXIS] = axis[X_AXIS]/abs(axis);
   axis[Y_AXIS] = axis[Y_AXIS]/abs(axis);
   axis[Z_AXIS] = axis[Z_AXIS]/abs(axis);

   double sina, cosa;
   if (angle == 90){
      sina = 1; cosa = 0;
   } else if (angle == 180){
      sina = 0; cosa = -1;
   } else {
      sina = sin(DEG_TO_RAD(angle));
      cosa = cos(DEG_TO_RAD(angle));
   }

   double uxsin = axis[X_AXIS]*sina;
   double uysin = axis[Y_AXIS]*sina;
   double uzsin = axis[Z_AXIS]*sina;
   
   double uxcos = (1-cosa)*axis[X_AXIS];
   double uycos = (1-cosa)*axis[Y_AXIS];
   double uzcos = (1-cosa)*axis[Z_AXIS];

   rotation[0][0] = axis[X_AXIS]*uxcos+cosa;
   rotation[0][1] = axis[Y_AXIS]*uxcos+uzsin;
   rotation[0][2] = axis[Z_AXIS]*uxcos-uysin;
   rotation[1][0] = axis[X_AXIS]*uycos-uzsin;
   rotation[1][1] = axis[Y_AXIS]*uycos+cosa;
   rotation[1][2] = axis[Z_AXIS]*uycos+uxsin;
   rotation[2][0] = axis[X_AXIS]*uzcos+uysin;
   rotation[2][1] = axis[Y_AXIS]*uzcos-uxsin;
   rotation[2][2] = axis[Z_AXIS]*uzcos+cosa;
}
