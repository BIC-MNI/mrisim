#ifndef __SAMPLE_H
#define __SAMPLE_H

/*****************************************************************************
 * 
 * SAMPLE.H
 *
 * Sample Class
 *
 * R. Kwan
 * August 3, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/sample.h,v 1.1 2003-05-30 16:43:14 bert Exp $
 * $Log: sample.h,v $
 * Revision 1.1  2003-05-30 16:43:14  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.2  1995/12/11  14:27:18  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.1  1995/08/08  14:45:36  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include <stdio.h>
#include "event.h"

/*****************************************************************************
 * Sample Class
 *****************************************************************************/

class Sample : public Event {
   public:
      Sample();                                // constructors
      Sample(Time_ms t);
      Sample(Time_ms t, void (*f)(Vector_3D&, void *), void *);
      Sample(const Sample& s);

      virtual Vector_3D& apply(Spin_Model& m);            // Sample the magnetization
      virtual void get_descriptor_string(char s[]);
      virtual Event *make_new_copy_of_event(void);

   protected:

      // Function to save samples
      void   (*_save_sample)(Vector_3D& m, void *);

      // Default save sample method
      static void _print_sample(Vector_3D& m, void *);

      // Data passed to save function 
      void   *_clientData;       
};

#endif

