/*****************************************************************************
 * 
 * TISSUE.CXX
 *
 * Tissue class member functions
 *
 * R. Kwan
 * March 15, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/tissue.cxx,v 1.1 2003-05-30 16:43:14 bert Exp $
// $Log: tissue.cxx,v $
// Revision 1.1  2003-05-30 16:43:14  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.2  1995/12/11  15:08:25  rkwan
// Updated for fast_iso_model.
//
// Revision 2.1  1995/08/08  14:55:36  rkwan
// Updated for new spin model and pulse sequence.
//

#include "tissue.h"

/*****************************************************************************
 * Tissue Class Methods
 *****************************************************************************/

/*****************************************************************************
 * Tissue Constructor
 *****************************************************************************/

Tissue::Tissue() {
   strcpy(_tissue_name,"NONE");
   _T1  = 0;
   _T2  = 0;
   _T2s = 0;
   _NH  = 0;
}

Tissue::Tissue(Time_ms T1, Time_ms T2, Time_ms T2s, float NH) {
   strcpy(_tissue_name,"NONE");
   _T1  = T1;
   _T2  = T2;
   _T2s = T2s;
   _NH  = NH;
}
 
Tissue::Tissue(const char *tissue_name, 
               Time_ms T1, Time_ms T2, Time_ms T2s, float NH) {
   strcpy(_tissue_name, tissue_name);
   _T1  = T1;
   _T2  = T2;
   _T2s = T2s;
   _NH  = NH;
}

Tissue::Tissue(const Tissue& t){
   strcpy(_tissue_name, t._tissue_name);
   _T1  = t._T1;
   _T2  = t._T2;
   _T2s = t._T2s;
   _NH  = t._NH;
}

/*****************************************************************************
 * Tissue::operator=
 * Assigns one tissue object to another.
 *****************************************************************************/

Tissue& Tissue::operator=(const Tissue& t){
   if (this != &t){
      strcpy(_tissue_name,t._tissue_name);
      _T1  = t._T1;
      _T2  = t._T2;
      _T2s = t._T2s;
      _NH  = t._NH;
   }
   return *this;
}

void Tissue::display_info(ostream& stream){
   stream << setw(20) << _tissue_name << " ";
   stream << setw(6) << _T1 << " "
          << setw(6) << _T2 << " "
          << setw(6) << _T2s << " "
          << setw(6) << _NH << " ";
}

void Tissue::dump_tissue_info(FILE *output){

   fprintf(output,"\n %s\n", _tissue_name);
   fprintf(output,"T1  = %7.1lf ms\n",_T1);
   fprintf(output,"T2  = %7.1lf ms\n",_T2);
   fprintf(output,"T2* = %7.1lf ms\n",_T2s);
   fprintf(output,"NH  = %4.2lf \n",_NH);

}



