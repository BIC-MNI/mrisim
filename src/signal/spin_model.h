#ifndef __SPIN_MODEL_H
#define __SPIN_MODEL_H

/*****************************************************************************
 * 
 * SPIN_MODEL.H
 *
 * Container class for MRI simulation spin models.
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/spin_model.h,v 1.2 2008-11-06 10:58:23 rotor Exp $
 * $Log: spin_model.h,v $
 * Revision 1.2  2008-11-06 10:58:23  rotor
 *  * fixed includes for iostream and friends
 *  * updated for new release (1.0.2)
 *
 * Revision 1.1  2003/05/30 16:43:14  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.4  1996/01/12  16:01:17  rkwan
 * Added flip error.
 *
 * Revision 2.3  1995/12/11  15:07:20  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.2  1995/08/08  18:46:46  rkwan
 * Added get_net_magnetization(float v[]) access function.
 *
 * Revision 2.1  1995/08/08  14:48:52  rkwan
 * MRI simulation spin model base class.
 *
 */

#include <iostream>
#include <mrisim/mrisim.h>
#include "tissue.h"
#include "vector.h"

/*****************************************************************************
 * Spin_Model Class
 *
 * Container class for MRI simulation spin models.
 * Defines the interface for simulation models, which should 
 * encapsulate information about 3D coordinates of spins and
 * tissue dependent characteristics.
 *
 *****************************************************************************/

class Spin_Model {
   public:
      Spin_Model();
       
      virtual ~Spin_Model();

      virtual void restore_equilibrium(void) = 0;
      virtual void rotate(Time_ms, Degrees, Axis) = 0;
      virtual void rotate(Time_ms, Degrees, Vector_3D&) = 0;
      virtual void relax(Time_ms) = 0;
      virtual void update(Time_ms) = 0;
      virtual void zero_transverse_magnetization(Time_ms) = 0;
      virtual void set_time(Time_ms) = 0;

      Vector_3D& get_net_magnetization(void);
      void       get_net_magnetization(float v[]);

      virtual Vector_3D& get_time_sample(Time_ms) = 0;

      operator Vector_3D&();

      double get_flip_error(void){
         return _flip_error; }
      void   set_flip_error(double flip_error){
         _flip_error = flip_error; }

   protected:
      Vector_3D _net_magnetization;

      // RF inhomogeneity flip angle error
      double _flip_error;

};

#endif

