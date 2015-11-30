#ifndef __MRISIMARGS_H
#define __MRISIMARGS_H

//==========================================================================
// MRISIMARGS.H
// Command line argument parsing for mrisim.
// Inherits from:  
// Base class to: 
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/mrisimargs.h,v 1.2 2004-08-10 15:39:07 bert Exp $
 * $Log: mrisimargs.h,v $
 * Revision 1.2  2004-08-10 15:39:07  bert
 * Add float.h
 *
 * Revision 1.1  2003/05/30 16:43:11  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.0.1.1  1996/07/16  20:42:24  rkwan
 * Command-line switch expansion.
 *
 * Revision 3.0  1996/06/18  19:57:33  rkwan
 * Fuzzy phantom update.
 *
 * Revision 2.5  1996/05/29  16:16:05  rkwan
 * Release 2.5
 *
 * Revision 1.3  1996/01/17  18:02:26  rkwan
 * Update for tx/rx coil modelling.
 *
 * Revision 1.2  1995/12/22  20:16:26  rkwan
 * Update for percent coil and RF map features.
 *
 * Revision 1.1  1995/12/11  15:25:01  rkwan
 * Initial revision
 *
 *========================================================================*/

#include "mrisim.h"
#include <limits.h>
#include <float.h>

extern "C" {
#include "ParseArgv.h"
};

//--------------------------------------------------------------------------
// Linked-list node for fuzzy phantom specifiers
//--------------------------------------------------------------------------

struct Fuzzy_Specifier {
   Fuzzy_Specifier() : filename(NULL), label(-1), next(NULL) {}
   ~Fuzzy_Specifier() {
      if (filename != NULL) free(filename);
   }

   char *filename; 
   int  label;
   Fuzzy_Specifier *next;
};

//--------------------------------------------------------------------------
// Phantom selection type
//--------------------------------------------------------------------------

enum Phantom_Type {DISCRETE, DISCRETE_RF, FUZZY, FUZZY_RF};

//--------------------------------------------------------------------------
// mrisimArgs class
// Command line argument parsing for mrisim.
//--------------------------------------------------------------------------

class mrisimArgs {
   public:
      mrisimArgs(int argc, char *argv[]);
      virtual ~mrisimArgs();

      // --- ParseArgv initialization table --- //

      static ArgvInfo argTable[];

      // --- Parameter and input files --- //

      static char *phantomFile;
      static char *tissueFile;
      static char *coilFile;
      static char *sequenceFile;
      static char *logFile;
      static char *rxmapFile;
      static char *txmapFile;
      static char *outputFile;

      static char *fuzzy_specifier_string;
      static Fuzzy_Specifier *fuzzy_list;  

      // --- Option flags --- //

      static int  logFlag;
      static int  clobberFlag;
      static int  verboseFlag;
      static int  versionFlag;
      static int  oldpvFlag;
      static int  oldoffsetFlag;

      // --- Sequence options --- //

      static double signal_gain;
      static double voxel_dims[3];
      static double voxel_offset[3];
      static int    recon_size;
      static int    matrix_size;
      static int    nslices;
      static double slice_thickness;

      // --- Access functions --- //

      inline int uses_fuzzy_phantom(void) const;
      inline int uses_discrete_phantom(void) const;
      inline int uses_default_fuzzy_phantom(void) const;
      inline int uses_default_discrete_phantom(void) const;
      inline int fuzzy_labels_specified(void) const;

      inline int uses_rx_map(void) const;
      inline int uses_tx_map(void) const;
      inline int uses_rf_phantom(void) const;
      inline int uses_default_rx_map(void) const;
      inline int uses_default_tx_map(void) const;

      // --- Convenience functions --- //

      Phantom_Type get_phantom_type(void) const;
      inline void delete_fuzzy_list(void);

   private:

      // --- Fuzzy phantom specifier list parsing --- //

      static int parse_fuzzy_list(char *head, int &count, Fuzzy_Specifier *list);
      static void delete_specifier_list(Fuzzy_Specifier *list);

};

//--------------------------------------------------------------------------
// ParseArgv helper function prototypes
//--------------------------------------------------------------------------

int  parse_optional_string_argument(char *dst, char *key, char *nextArg);

//--------------------------------------------------------------------------
// Inline member functions
//--------------------------------------------------------------------------

//--------------------------------------------------------------------------
// mrisimArgs::uses_fuzzy_phantom
// Returns TRUE if -fuzzy was used.
//--------------------------------------------------------------------------

inline 
int mrisimArgs::uses_fuzzy_phantom(void) const {
   return (mrisimArgs::fuzzy_specifier_string != NULL);
}

//--------------------------------------------------------------------------
// mrisimArgs::uses_discrete_phantom
// Returns TRUE if -discrete was used.
//--------------------------------------------------------------------------

inline 
int mrisimArgs::uses_discrete_phantom(void) const {
   return (mrisimArgs::phantomFile != NULL);
}

//--------------------------------------------------------------------------
// mrisimArgs::uses_default_fuzzy_phantom
// Returns TRUE if -fuzzy was used without specifying phantom files.
//--------------------------------------------------------------------------

inline 
int mrisimArgs::uses_default_fuzzy_phantom(void) const {
   return (uses_fuzzy_phantom() && 
           (mrisimArgs::fuzzy_specifier_string[0]=='\0'));
}

//--------------------------------------------------------------------------
// mrisimArgs::uses_default_discrete_phantom
// Returns TRUE if -discrete was used without specifying a phantom file.
//--------------------------------------------------------------------------

inline
int mrisimArgs::uses_default_discrete_phantom(void) const {
   return (uses_discrete_phantom() &&
           (mrisimArgs::phantomFile[0]=='\0'));
}

//--------------------------------------------------------------------------
// mrisimArgs::fuzzy_labels_specified
// Returns TRUE if -fuzzy was used with labels specified in the phantom
// list.
//--------------------------------------------------------------------------

inline
int mrisimArgs::fuzzy_labels_specified(void) const {
   return (uses_fuzzy_phantom() && 
           !uses_default_fuzzy_phantom() && 
           (fuzzy_list->label != -1));
} 

//--------------------------------------------------------------------------
// mrisimArgs::uses_rx_map
// Returns TRUE if -rxmap was used.
//--------------------------------------------------------------------------

inline
int mrisimArgs::uses_rx_map(void) const {
   return (mrisimArgs::rxmapFile != NULL);
}

//--------------------------------------------------------------------------
// mrisimArgs::uses_tx_map
// Returns TRUE if -txmap was used.
//--------------------------------------------------------------------------

inline
int mrisimArgs::uses_tx_map(void) const {
   return (mrisimArgs::txmapFile != NULL);
}

//--------------------------------------------------------------------------
// mrisimArgs::uses_default_rx_map
// Returns TRUE if -rxmap was used without specifying a map file.
//--------------------------------------------------------------------------

inline
int mrisimArgs::uses_default_rx_map(void) const {
   return (uses_rx_map() &&
           (mrisimArgs::rxmapFile[0]=='\0'));
}

//--------------------------------------------------------------------------
// mrisimArgs::uses_default_tx_map
// Returns TRUE if -txmap was used without specifying a map file.
//--------------------------------------------------------------------------

inline
int mrisimArgs::uses_default_tx_map(void) const {
   return (uses_tx_map() &&
           (mrisimArgs::txmapFile[0]=='\0'));
}

//--------------------------------------------------------------------------
// mrisimArgs::uses_rf_phantom
// Returns TRUE if RF inhomogeneity is to be simulated.
//--------------------------------------------------------------------------

inline
int mrisimArgs::uses_rf_phantom(void) const {
   return (uses_rx_map() || uses_tx_map());
}

//--------------------------------------------------------------------------
// mrisimArgs::delete_fuzzy_list
// Deletes the linked-list of fuzzy specifiers.
//--------------------------------------------------------------------------

inline
void mrisimArgs::delete_fuzzy_list(void) {
   if (mrisimArgs::fuzzy_list != NULL) {
      delete_specifier_list(mrisimArgs::fuzzy_list);
      mrisimArgs::fuzzy_list = NULL;
   }
}

#endif
 
