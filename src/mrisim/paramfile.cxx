//==========================================================================
// PARAMFILE.CXX
// Parameter file parsing base class for mrisim.
// Inherits from:  ifstream
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/paramfile.cxx,v 1.1 2003-05-30 16:43:11 bert Exp $
// $Log: paramfile.cxx,v $
// Revision 1.1  2003-05-30 16:43:11  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  16:05:17  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:17:05  rkwan
// Release 2.5
//
// Revision 1.2  1995/12/11  15:29:38  rkwan
// Doc update.
//
//==========================================================================

#include "paramfile.h"
#include <stdlib.h>

//--------------------------------------------------------------------------
// ParamFile::_tmp
// static string storage for parameter file input.
//--------------------------------------------------------------------------

const unsigned int ParamFile::_tmp_length = 255;
char ParamFile::_tmp[255];

//--------------------------------------------------------------------------
// ParamFile constructor
//--------------------------------------------------------------------------

ParamFile::ParamFile(const char *path) : ifstream(path) {}

//--------------------------------------------------------------------------
// ParamFile destructor
//--------------------------------------------------------------------------

ParamFile::~ParamFile() {
   this->close();
}

//--------------------------------------------------------------------------
// ParamFile::ignoreline
// Ignores a line of input in the parameter file up to the first 
// occurence of end_marker.  By default, end_marker is a newline.
//--------------------------------------------------------------------------

istream& ParamFile::ignoreline(char end_marker){
   this->ignore(132,end_marker);
   return *this;
}

//--------------------------------------------------------------------------
// ParamFile::getfield
// Reads a field from the parameter file.  
// End of fields are marked by newlines.
// Overloaded for the required data type.
//--------------------------------------------------------------------------

istream& ParamFile::getfield(unsigned char &uchar){
   this->getfield(_tmp_length, _tmp);
   uchar = (unsigned char)atoi(_tmp);
   return *this;
}

istream& ParamFile::getfield(int& i){
   this->getfield(_tmp_length, _tmp);
   i = atoi(_tmp);
   return *this;
}

istream& ParamFile::getfield(long& l){
   this->getfield(_tmp_length, _tmp);
   l = atol(_tmp);
   return *this;
}

istream& ParamFile::getfield(double& d){
   this->getfield(_tmp_length, _tmp);
   d = atof(_tmp);
   return *this;
}

istream& ParamFile::getfield(unsigned int length, char str[]){

   enum State {LABEL=0, VALUE=1, COMMENT=2, END_COMMENT=3};

   State state = LABEL;
   char  c, *cptr = str;

   while(this->get(c) && cptr < str+length) {
      switch(state) {
         case LABEL:
            if (c == '#') {
               state = COMMENT;
            } else if (c == ':') {
               state = VALUE;
            }
            break; 
         case VALUE:
            if (c == '#' || c == '(') {
               state = END_COMMENT;
               *cptr++ = '\0';
            } else if (c == '\n') {
               *cptr++ = '\0';
               return *this;
            } else if (c != '\t' && c != ' ') {
               *cptr++ = c;
            }
            break;
         case COMMENT:
            if (c == '\n') {
               state = LABEL;
            } 
            break;
         case END_COMMENT:
            if (c == '\n') {
               return *this;
            } 
      }
   }
   *cptr++ = '\0';
   return *this;
}
 
 
istream& ParamFile::getlist2(double &param1, double &param2) {

   getfield(255, _tmp);

   char *cptr = _tmp;
   char value[255], *vptr = value; 

   while(*cptr != ',' && *cptr != '\0') {
      *vptr++ = *cptr++;
   }
   if (*cptr == ',') cptr++;
   *vptr++ = '\0';
   param1 = atof(value);

   vptr = value;
   while(*cptr != '\0') {
      *vptr++ = *cptr++;
   }
   *vptr++ = '\0';
   param2 = atof(value);

   return *this;
}
