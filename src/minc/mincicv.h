#ifndef __MINCICV_H
#define __MINCICV_H

//==========================================================================
// MINCICV.H
// MINC ICV class MINC_ICV.
// Inherits from:  
// Base class to:  MINC_File
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/mincicv.h,v 1.2 2008-11-06 10:58:23 rotor Exp $
 * $Log: mincicv.h,v $
 * Revision 1.2  2008-11-06 10:58:23  rotor
 *  * fixed includes for iostream and friends
 *  * updated for new release (1.0.2)
 *
 * Revision 1.1  2003/05/30 16:43:09  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:27:44  rkwan
 * Release 2.5
 *
 * Revision 1.4  1995/12/15  21:07:11  rkwan
 * Added set_default_float_icv.
 *
 * Revision 1.3  1995/12/11  15:10:11  rkwan
 * Fix RCS header bug.
 *
 *=========================================================================*/

#include <iostream>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

extern "C" {
#include <minc.h>
#undef public
#define TRUE 1
#define FALSE 0
}

#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#endif

#ifndef NC_OPTS_VAL
#define NC_OPTS_VAL NC_VERBOSE | NC_FATAL
#endif

//--------------------------------------------------------------------------
// MINC_ICV class
// Describes a MINC icv which controls the appearance of MINC volumes.
//--------------------------------------------------------------------------

class MINC_ICV {
   public:
      MINC_ICV();

      virtual ~MINC_ICV();

      // MINC ICV property access functions

      inline int set_icv_property(int icv_property, double value);
      inline int set_icv_property(int icv_property, int value);
      inline int set_icv_property(int icv_property, long value);
      inline int set_icv_property(int icv_property, char *value);
 
      inline int get_icv_property(int icv_property, double *value);
      inline int get_icv_property(int icv_property, int *value);
      inline int get_icv_property(int icv_property, long *value);
      inline int get_icv_property(int icv_property, char *value);

      inline int get_icvid(void) const {return _icvid;}
   
      // Convenience member functions

      const char *get_icv_type_name(void);
      void set_default_byte_icv(void); 
      void set_default_short_icv(void);
      void set_default_float_icv(void);

  
   protected:
      
      int _icvid;

};
 
//--------------------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// MINC_ICV::set_icv_property
// Sets named icv properties.  Overloaded for the required data type.
//--------------------------------------------------------------------------

inline
int MINC_ICV::set_icv_property(int icv_property, double value) {
   return miicv_setdbl(_icvid, icv_property, value);
}

inline
int MINC_ICV::set_icv_property(int icv_property, int value) {
   return miicv_setint(_icvid, icv_property, value);
}

inline
int MINC_ICV::set_icv_property(int icv_property, long value) {
   return miicv_setlong(_icvid, icv_property, value);
}

inline
int MINC_ICV::set_icv_property(int icv_property, char *value) {
   return miicv_setstr(_icvid, icv_property, value);
}

//--------------------------------------------------------------------------
// MINC_ICV::get_icv_property
// Inquires about named icv properties.  Overloaded for the required
// data type.
//--------------------------------------------------------------------------

inline
int MINC_ICV::get_icv_property(int icv_property, double *value){

#ifdef DEBUG
   assert(value != NULL);
#endif
   return miicv_inqdbl(_icvid, icv_property, value);
}

inline
int MINC_ICV::get_icv_property(int icv_property, int *value){

#ifdef DEBUG
   assert(value != NULL);
#endif

   return miicv_inqint(_icvid, icv_property, value);
}

inline
int MINC_ICV::get_icv_property(int icv_property, long *value){

#ifdef DEBUG
   assert(value != NULL);
#endif

   return miicv_inqlong(_icvid, icv_property, value);
}

inline
int MINC_ICV::get_icv_property(int icv_property, char *value){

#ifdef DEBUG
   assert(value != NULL);
#endif

   return miicv_inqstr(_icvid, icv_property, value);
}



#endif
