//===========================================================================
// MRIMATRIX.CXX
//
// R.Kwan
// August 31, 1995
//
// (C) Copyright 1995 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/mristring.cxx,v 1.1 2003-05-30 16:43:09 bert Exp $
// $Log: mristring.cxx,v $
// Revision 1.1  2003-05-30 16:43:09  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:28:30  rkwan
// Release 2.5
//
// Revision 1.2  1995/12/11  14:17:52  rkwan
// Updated for fast_iso_model.
//
//===========================================================================

#include "mristring.h"

//--------------------------------------------------------------------------
// MRI_String constructors
//--------------------------------------------------------------------------

MRI_String::MRI_String() {
   _length = DEFAULT_STRING_LENGTH;
   _string = new char[_length+1];
}

MRI_String::MRI_String(unsigned int length) {
   _length = length;
   _string = new char[_length+1];
}

MRI_String::MRI_String(const char *string) {
   _length = strlen(string);
   _string = new char[_length+1];
   strcpy(_string, string);
}

MRI_String::MRI_String(const MRI_String& string) {
   if (this != &string){
      delete[] _string;
      _length = string._length;
      _string = new char[_length+1];
      strcpy(_string, string._string);
   }
}

//--------------------------------------------------------------------------
// MRI_String destructor
//--------------------------------------------------------------------------

MRI_String::~MRI_String() {
   delete[] _string;
}

//--------------------------------------------------------------------------
// MRI_String::operator=
// Assignment operator
//--------------------------------------------------------------------------

MRI_String& MRI_String::operator=(const char *string) {
   int length;
   if ((length = strlen(string)) > _length){
      delete[] _string;
      _length = length;
      _string = new char[_length+1];
   }
   strcpy(_string, string);
   return *this;
}

MRI_String& MRI_String::operator=(const MRI_String& string) {
   if (this != &string){
      if (string._length > _length){
         delete[] _string;
         _length = string._length;
         _string = new char[_length+1];
      }
      strcpy(_string, string._string);
   }
   return *this;

}

//--------------------------------------------------------------------------
// MRI_String::operator()
// Returns a substring of the string.
//--------------------------------------------------------------------------

MRI_String MRI_String::operator() (unsigned int n1, unsigned int n2) {
   MRI_String tmp(n2-n1+1);
   strncpy(tmp._string, &(_string[n1]), n2-n1+1);
   return tmp;
}

