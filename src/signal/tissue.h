#ifndef __TISSUE_H
#define __TISSUE_H

/*****************************************************************************
 * 
 * TISSUE.H
 *
 * Tissue classes
 *
 * R. Kwan
 * March 15, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/tissue.h,v 1.2 2008-11-06 10:58:23 rotor Exp $
 * $Log: tissue.h,v $
 * Revision 1.2  2008-11-06 10:58:23  rotor
 *  * fixed includes for iostream and friends
 *  * updated for new release (1.0.2)
 *
 * Revision 1.1  2003/05/30 16:43:14  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.3  1995/12/11  15:08:17  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.1  1995/08/08  14:55:30  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include <mrisim/mrisim.h>
#include <minc/mristring.h>
#include <iostream>
#include <iomanip>

using namespace std;

/*****************************************************************************
 * Tissue Class
 *
 * Encapsulates intrinsic tissue parameters
 *   T1  longitudinal relaxation time (ms)
 *   T2  transverse relaxation time (ms)
 *   T2s T2* transverse relaxation time (ms)
 *   NH  relative proton density
 *
 * Can be extended to include susceptibility, chemical shift, etc...
 *
 *****************************************************************************/

class Tissue {
   public:
      Tissue();                                // constructors
      Tissue(Time_ms T1, Time_ms T2, 
             Time_ms T2s, float NH);
      Tissue(const char *tissue_name, 
             Time_ms T1, Time_ms T2, Time_ms T2s, float NH);
      Tissue(const Tissue& t);

      Tissue& operator=(const Tissue& t);      // assignment

      Time_ms get_T1(void) const {return _T1;} // access functions
      Time_ms get_T2(void) const {return _T2;}
      Time_ms get_T2s(void) const {return _T2s;}
      float   get_NH(void) const {return _NH;}
      char *get_tissue_name(void) {return _tissue_name;}

      void  set_T1(Time_ms T1) {_T1 = T1;}
      void  set_T2(Time_ms T2) {_T2 = T2;}
      void  set_T2s(Time_ms T2s) {_T2s = T2s;}
      void  set_NH(float  NH) {_NH = NH;}

      void  display_info(ostream&  stream);
      void  dump_tissue_info(FILE *output);

   protected:
      char  _tissue_name[80];                  // tissue description
      Time_ms   _T1;                           // longitudinal relaxation (ms)
      Time_ms   _T2;                           // transverse relaxation (ms)
      Time_ms   _T2s;                          // T2* trans. relaxation (ms)
      float     _NH;                           // relative proton density
};

#endif

