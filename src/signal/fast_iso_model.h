#ifndef __FAST_ISO_MODEL_H
#define __FAST_ISO_MODEL_H

/*****************************************************************************
 * 
 * FAST_ISO_MODEL.H
 *
 * Fast Isochromat spin system model uses FFT.
 *
 * R. Kwan
 * November 22, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/fast_iso_model.h,v 1.1 2003-05-30 16:43:12 bert Exp $
 * $Log: fast_iso_model.h,v $
 * Revision 1.1  2003-05-30 16:43:12  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:31:50  rkwan
 * Release 2.5
 *
 * Revision 1.1  1995/12/11  14:23:25  rkwan
 * Initial revision
 *
 */

#include <math.h>
#ifdef DEBUG
#include <assert.h>
#endif

#include <mrisim/mrisim.h>
#include "tissue.h"
#include "spin_model.h"

extern "C" {
#include "../minc/fourn.h"
}

/*****************************************************************************
 * Fast_Isochromat_Model Class
 *****************************************************************************/

class Fast_Isochromat_Model : public Spin_Model, public Tissue {
   public:
      Fast_Isochromat_Model(int n);
      Fast_Isochromat_Model(Time_ms T1, Time_ms T2, Time_ms T2s, float NH, 
                            int n=0);
      Fast_Isochromat_Model(const Tissue& tissue, int n=0);

      virtual ~Fast_Isochromat_Model();

      virtual void restore_equilibrium(void);
      virtual void rotate(Time_ms t, Degrees angle, Axis axis);
      virtual void rotate(Time_ms t, Degrees angle, Vector_3D& axis);
      virtual void relax(Time_ms t);
      virtual void update(Time_ms t);
      virtual void set_time(Time_ms t);
      virtual void zero_transverse_magnetization(Time_ms t);
      virtual void spoil(double G, Time_ms t);

      void use_linear_resonant_freq(Hertz bandwidth);
      void use_lorentzian_distribution(double alpha);
      void use_uniform_distribution(void);
      void use_exponential_decay_model(void);
      void use_custom_resonant_freq(double freq[]);
      void use_custom_distribution(double dist[]);

      int     get_num_of_isochromats(void){
         return _num_of_isochromats;}
      Hertz   get_bandwidth(void){
         return _bandwidth;}
      Time_ms get_time_step(void) {
         return (Time_ms)(1.0/(2*_bandwidth)); }
      void    get_magnetization(int n, float v[]);

      void get_time_samples(double xy_samples[], double z_samples[]);
      void get_time_samples(double xy_samples[], double z_samples[],
                            int num_of_samples, Time_ms t_start,
                            Time_ms t_step);
      void get_time_sample(Vector_3D& sample, int n);
      Vector_3D& get_time_sample(Time_ms t);

      void relax_slow(Time_ms t);

      virtual void display_model_info(ostream& stream);
      virtual void display_model_state(ostream& stream);

   protected:
      int       _num_of_isochromats;
      Time_ms   _t_0;
      Time_ms   _t;
 
   private:
      double    *_m_equil;
      Hertz     *_off_resonance_freq;
      double    *_m_0_z;
      double    *_m_0_xy;
      double    *_m_z;
      double    *_m_xy;
      Hertz     _bandwidth;

      void _compute_net_mag(Vector_3D& net, double z_samples[], 
                            double xy_samples[]);
      void _update_initial_mag(Time_ms t);
      void _compute_rotation_matrix(double rotation[][3],
                                    Axis axis, Degrees angle);
      void _compute_general_rotation_matrix(double rotation[][3],
                                    Vector_3D& axis, Degrees angle);
};

#endif

