#ifndef __VECTOR_MODEL_H
#define __VECTOR_MODEL_H

/*****************************************************************************
 * 
 * VECTOR_MODEL.H
 *
 * Single magnetization vector model.
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/vector_model.h,v 1.1 2003-05-30 16:43:14 bert Exp $
 * $Log: vector_model.h,v $
 * Revision 1.1  2003-05-30 16:43:14  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 1.1  1995/12/11  14:32:04  rkwan
 * Initial revision
 *
 */

#include <math.h>
#ifdef DEBUG
#include <assert.h>
#endif

#include "tissue.h"
#include "spin_model.h"

/*****************************************************************************
 * Vector_Model Class
 *****************************************************************************/

class Vector_Model : public Spin_Model, public Tissue {
   public:
      Vector_Model();
      Vector_Model(Time_ms T1, Time_ms T2, Time_ms T2s, float NH);
      Vector_Model(const Tissue& tissue);

      virtual ~Vector_Model() {}

      virtual void restore_equilibrium(void);
      virtual void rotate(Time_ms t, Degrees angle, Axis axis);
      virtual void rotate(Time_ms t, Degrees angle, Vector_3D& axis);
      virtual void relax(Time_ms t);
      virtual void update(Time_ms t);
      virtual void set_time(Time_ms t);
      virtual void zero_transverse_magnetization(Time_ms t);

      virtual void display_model_info(ostream& stream);
      virtual void display_model_state(ostream& stream);
   
      Vector_3D& get_time_sample(Time_ms t){
         this->relax(t);
         return (Vector_3D&)*this;
      }
   
   protected:
      double    _m_equil;
      Vector_3D _m_0;
      Time_ms   _t_0;
      Time_ms   _t;

};

#endif

