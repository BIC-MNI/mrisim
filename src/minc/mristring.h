#ifndef __MRISTRING_H
#define __MRISTRING_H

//===========================================================================
// MRISTRING.H
//
// R.Kwan
// August 31, 1995
//
// (C) Copyright 1995 by R.Kwan
//===========================================================================

/*===========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/mristring.h,v 1.1 2003-05-30 16:43:09 bert Exp $
 * $Log: mristring.h,v $
 * Revision 1.1  2003-05-30 16:43:09  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:28:23  rkwan
 * Release 2.5
 *
 * Revision 1.3  1995/12/11  14:36:49  rkwan
 * Fix RCS header.
 *
 *=========================================================================*/

#include <string.h>

#ifdef DEBUG
#include <assert.h>
#endif

#ifndef DEFAULT_STRING_LENGTH
#define DEFAULT_STRING_LENGTH 20
#endif

//---------------------------------------------------------------------------
// MRI_String class
//---------------------------------------------------------------------------

class MRI_String {
   public:
      MRI_String();
      MRI_String(unsigned int length);
      MRI_String(const char *string);
      MRI_String(const MRI_String& string);

      virtual ~MRI_String();

      MRI_String& operator=(const char *string);
      MRI_String& operator=(const MRI_String& string);

      unsigned int get_string_length(void) {
         return _length; }

      char& operator[] (unsigned int n){
         return _string[n]; }
      operator char *() {
         return _string; }
      MRI_String operator() (unsigned int n1, unsigned int n2);

   protected:
      char *_string;
      int  _length;

};

#endif


   
