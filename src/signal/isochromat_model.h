#ifndef __ISOCHROMAT_MODEL_H
#define __ISOCHROMAT_MODEL_H

/*****************************************************************************
 * 
 * ISOCHROMAT_MODEL.H
 *
 * Isochromat spin system model class.
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/isochromat_model.h,v 1.1 2003-05-30 16:43:13 bert Exp $
 * $Log: isochromat_model.h,v $
 * Revision 1.1  2003-05-30 16:43:13  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:32:58  rkwan
 * Release 2.5
 *
 * Revision 2.4  1996/01/17  17:58:17  rkwan
 * Consistency fix between Fast_Isochromat_Model and Isochromat_Model.
 *
 * Revision 2.3  1995/12/11  15:08:34  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.2  1995/12/11  14:23:54  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.1  1995/08/08  14:37:54  rkwan
 * Updated for new spin model and pulse sequence.
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
 * Isochromat_Model Class
 *****************************************************************************/

class Isochromat_Model : public Spin_Model, public Tissue {
   public:
      Isochromat_Model(int n);
      Isochromat_Model(Time_ms T1, Time_ms T2, Time_ms T2s, float NH, 
                            int n=0);
      Isochromat_Model(const Tissue& tissue, int n=0);

      virtual ~Isochromat_Model();

      virtual void restore_equilibrium(void);
      virtual void rotate(Time_ms t, Degrees angle, Axis axis);
      virtual void rotate(Time_ms t, Degrees angle, Vector_3D& axis);
      virtual void relax(Time_ms t);
      virtual void update(Time_ms t);
      virtual void set_time(Time_ms t);
      virtual void zero_transverse_magnetization(Time_ms t);

      void use_linear_resonant_freq(Hertz bandwidth);
      void use_lorentzian_distribution(double alpha);
      void use_uniform_distribution(void);
      void use_exponential_decay_model(void);
      void use_custom_resonant_freq(double freq[]);
      void use_custom_distribution(double dist[]);

      int     get_num_of_isochromats(void){
         return _num_of_isochromats; }
      Hertz   get_bandwidth(void){ 
         return _bandwidth; }
      Time_ms get_time_step(void){
         return (Time_ms)(1.0/(2*_bandwidth)); }
      void   get_magnetization(int n, float v[]);

      Vector_3D& get_time_sample(Time_ms t){
         this->relax(t);
         return (Vector_3D&)*this;
      }

      virtual void display_model_info(ostream& stream);
      virtual void display_model_state(ostream& stream);

   protected:
      int       _num_of_isochromats;
      Time_ms   _t_0;
      Time_ms   _t;

   private:
      double    *_m_equil;
      Hertz     *_off_resonance_freq;
      Vector_3D *_m_0;
      Vector_3D *_m;
      Hertz     _bandwidth;

};

#endif

