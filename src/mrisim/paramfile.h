#ifndef __PARAMFILE_H
#define __PARAMFILE_H

//==========================================================================
// PARAMFILE.H
// Parameter file parsing base class for mrisim.
// Inherits from:  ifstream
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/paramfile.h,v 1.3 2008-11-06 10:58:23 rotor Exp $
 * $Log: paramfile.h,v $
 * Revision 1.3  2008-11-06 10:58:23  rotor
 *  * fixed includes for iostream and friends
 *  * updated for new release (1.0.2)
 *
 * Revision 1.2  2004/08/10 15:38:05  bert
 * Add iostream.h
 *
 * Revision 1.1  2003/05/30 16:43:11  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  16:05:02  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:16:56  rkwan
 * Release 2.5
 *
 * Revision 1.2  1995/12/11  15:29:26  rkwan
 * Doc update.
 *
 *========================================================================*/
#include <iostream>
#include <fstream>

using namespace std;

//--------------------------------------------------------------------------
// ParamFile class
// Parameter file parsing base class for mrisim.
//--------------------------------------------------------------------------

class ParamFile : public ifstream {
   public:
      ParamFile(const char *path);

      virtual ~ParamFile();

      istream& ignoreline(char end_marker = '\n');
      istream& getfield(unsigned char& uchar); 
      istream& getfield(int& i);
      istream& getfield(long& l);
      istream& getfield(double& d);
      istream& getfield(unsigned int length, char str[]);

      istream& getlist2(double &param1, double &param2);

   protected:
      static const unsigned int _tmp_length;
      static char _tmp[255];
};

#endif
