#ifndef __QUICK_MODEL_H
#define __QUICK_MODEL_H

/*****************************************************************************
 * 
 * QUICK_MODEL.H
 *
 * Quick spin model for use with precalculated pulse sequences.
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/quick_model.h,v 1.1 2003-05-30 16:43:13 bert Exp $
 * $Log: quick_model.h,v $
 * Revision 1.1  2003-05-30 16:43:13  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.4  1995/12/11  15:51:11  rkwan
 * Added arbitrary axis rotate.
 *
 * Revision 2.3  1995/12/11  15:46:47  rkwan
 * Add get_time_sample.
 *
 * Revision 2.2  1995/12/11  15:07:01  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.1  1995/08/08  14:41:15  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include "tissue.h"
#include "spin_model.h"

/*****************************************************************************
 * Quick_Model Class
 *****************************************************************************/

class Quick_Model : public Spin_Model, public Tissue {
   public:
      Quick_Model();
      Quick_Model(Time_ms T1, Time_ms T2, Time_ms T2s, float NH);
      Quick_Model(const Tissue& tissue);

      ~Quick_Model() { }

      void restore_equilibrium(void) { }
      void rotate(Time_ms, Degrees, Axis) { }
      void rotate(Time_ms, Degrees, Vector_3D&) { }
      void relax(Time_ms) { }
      void update(Time_ms) { }
      void zero_transverse_magnetization(Time_ms) { }
      void set_time(Time_ms) { }
      
      Vector_3D& get_time_sample(Time_ms) {return _net_magnetization;}

      void set_net_magnetization(Vector_3D& v); 
};

#endif

