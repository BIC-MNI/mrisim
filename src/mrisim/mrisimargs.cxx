//==========================================================================
// MRISIMARGS.CXX
// Command line argument parsing for mrisim.
// Inherits from:  
// Base class to: 
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/mrisimargs.cxx,v 1.1 2003-05-30 16:43:11 bert Exp $
// $Log: mrisimargs.cxx,v $
// Revision 1.1  2003-05-30 16:43:11  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  19:28:44  rkwan
// Release 3.1 update.
//
// Revision 3.0.2.1  1996/07/17  15:08:57  rkwan
// Maintenance release:  extended command-line switches.
//
// Revision 3.0.1.1  1996/07/08  14:32:25  rkwan
// Maintenance release: Fuzzy normalization performed by phantom.
//
// Revision 3.0  1996/06/18  20:29:45  rkwan
// Fuzzy phantom update.
//
// Revision 2.5  1996/05/29  16:16:10  rkwan
// Release 2.5
//
// Revision 1.3  1996/01/17  18:02:18  rkwan
// Update for tx/rx coil modelling.
//
// Revision 1.2  1995/12/22  20:16:33  rkwan
// Update for percent coil and RF map features.
//
// Revision 1.1  1995/12/11  15:25:18  rkwan
// Initial revision
//
//==========================================================================

#include <string.h>
#include <ctype.h>

#include "mrisimargs.h"

//------------------------------------------------------------------------- 
// Static command line argument flags
//------------------------------------------------------------------------- 

// --- Parameter and input files --- //

char *mrisimArgs::phantomFile  = NULL;
char *mrisimArgs::tissueFile   = NULL;
char *mrisimArgs::coilFile     = NULL;
char *mrisimArgs::sequenceFile = NULL;
char *mrisimArgs::logFile      = NULL;
char *mrisimArgs::rxmapFile    = NULL;
char *mrisimArgs::txmapFile    = NULL;
char *mrisimArgs::outputFile   = NULL;

char            *mrisimArgs::fuzzy_specifier_string = NULL;
Fuzzy_Specifier *mrisimArgs::fuzzy_list             = NULL;

// --- Option flags --- //

int mrisimArgs::logFlag        = FALSE;
int mrisimArgs::clobberFlag    = FALSE;
int mrisimArgs::verboseFlag    = TRUE;
int mrisimArgs::versionFlag    = FALSE;
int mrisimArgs::oldpvFlag      = FALSE;
int mrisimArgs::oldoffsetFlag  = FALSE;

// --- Sequence options --- //

double mrisimArgs::signal_gain     = 10000;
double mrisimArgs::voxel_dims[3]   = {0,0,0};
double mrisimArgs::voxel_offset[3] = {DBL_MIN, DBL_MIN, DBL_MIN};
int    mrisimArgs::recon_size      = 0;
int    mrisimArgs::matrix_size     = 0;
int    mrisimArgs::nslices         = 0;
double mrisimArgs::slice_thickness = 0;

//------------------------------------------------------------------------- 
// Command line argument descriptor table
//------------------------------------------------------------------------- 

ArgvInfo mrisimArgs::argTable[] = {
   {NULL, ARGV_HELP, (char *)NULL, (char *) NULL,
            ""},
   {"-noclobber", ARGV_CONSTANT, (char *)FALSE, 
            (char *)&mrisimArgs::clobberFlag,
            "Do not overwrite output file (default)."},
   {"-clobber", ARGV_CONSTANT, (char *)TRUE, 
            (char *)&mrisimArgs::clobberFlag,
            "Overwrite output file."},
   {"-verbose", ARGV_CONSTANT, (char *)TRUE, 
            (char *)&mrisimArgs::verboseFlag,
            "Print out log messages (default)."},
   {"-quiet", ARGV_CONSTANT, (char *)FALSE, 
            (char *)&mrisimArgs::verboseFlag,
            "Do not print out log messages."},
   {"-version", ARGV_CONSTANT, (char *)TRUE,
            (char *)&mrisimArgs::versionFlag,
            "Print version information and exit."},
   {"-tissue", ARGV_STRING, (char *)1, 
            (char *)&mrisimArgs::tissueFile,
            "Tissue parameter file."},
   {"-coil", ARGV_STRING, (char *)1, 
            (char *)&mrisimArgs::coilFile,
            "RF Coil parameter file."},
   {"-sequence", ARGV_STRING, (char *)1, 
            (char *)&mrisimArgs::sequenceFile,
            "Pulse sequence parameter file."},
   {"-discrete", ARGV_FUNC, 
            (char *)parse_optional_string_argument, 
            (char *)&mrisimArgs::phantomFile,
            "Use a discrete phantom."},
   {"-fuzzy", ARGV_FUNC, 
            (char *)parse_optional_string_argument,
            (char *)&mrisimArgs::fuzzy_specifier_string,
            "Use a fuzzy phantom."},
   {"-logfile", ARGV_FUNC, 
            (char *)parse_optional_string_argument,
            (char *)&mrisimArgs::logFile,
            "Name of log file to generate."},
   {"-rxmap", ARGV_FUNC,
            (char *)parse_optional_string_argument,
            (char *)&mrisimArgs::rxmapFile,
            "Model receive coil inhomogeneity."},
   {"-txmap", ARGV_FUNC,
            (char *)parse_optional_string_argument,
            (char *)&mrisimArgs::txmapFile,
            "Model transmit coil inhomogeneity."},
   {"-nnpv", ARGV_CONSTANT, (char *)TRUE, 
            (char *)&mrisimArgs::oldpvFlag,
            "Use nearest neighbour partial volume evaluation."},
   {"-old_offset", ARGV_CONSTANT, (char *)TRUE,
            (char *)&mrisimArgs::oldoffsetFlag,
            "Use old quarter-voxel shift offset when resampling."},
   {"-gain", ARGV_FLOAT, (char *) 1, 
            (char *)&mrisimArgs::signal_gain,
            "Output image signal gain multiplier."},
   {"-voxel_dims", ARGV_FLOAT, (char *) 3, 
            (char *)&mrisimArgs::voxel_dims,
            "Voxel dimensions of simulated image (mm) in (z y x) order."},
   {"-voxel_offset", ARGV_FLOAT, (char *) 3,
            (char *)&mrisimArgs::voxel_offset,
            "Voxel offset of the simulated image (mm) in (z y x) order."},
   {"-recon_matrix", ARGV_INT, (char *) 1,
            (char *)&mrisimArgs::recon_size,
            "Specify the reconstructed image matrix size."},
   {"-scan_matrix", ARGV_INT, (char *) 1,
            (char *)&mrisimArgs::matrix_size,
            "Specify the acquisition matrix size."},
   {"-slice_thickness", ARGV_FLOAT, (char *) 1,
            (char *)&mrisimArgs::slice_thickness,
            "Specify the acquisition slice thickness (mm)."},
   {"-nslices", ARGV_INT, (char *) 1,
             (char *)&mrisimArgs::nslices,
             "Specify the acquisition slice thickness (mm)."},
   {(char *)NULL, ARGV_END, (char *)NULL, (char *)NULL,
            (char *)NULL}
};

//------------------------------------------------------------------------- 
// mrisimArgs constructor
//------------------------------------------------------------------------- 

mrisimArgs::mrisimArgs(int argc, char *argv[]){

   // Parse and check command line arguments

   if (ParseArgv(&argc, argv, argTable, 0) || (argc != 2) ||
       mrisimArgs::versionFlag){

      if (mrisimArgs::versionFlag){
         cerr << endl << argv[0] << " $State: Exp $ $Revision: 1.1 $" 
              << endl;
         cerr << "Build Date: $Date: 2003-05-30 16:43:11 $" << endl << endl;
         exit(0);
      } else {
         cerr << endl << "Usage: " << argv[0] << " [<options>] <output.mnc>" 
              << endl;
         cerr <<         "       " << argv[0] << " [-help]" << endl << endl;
         exit(EXIT_FAILURE);
      }

   }
   
   // Check output file
   mrisimArgs::outputFile = argv[1];
   if (strcmp(&mrisimArgs::outputFile[strlen(argv[1])-4],".mnc") != 0){
      cerr << "Output file name must end in extension .mnc" << endl;
      exit(EXIT_FAILURE);
   }

   // Ensure all required parameter files are given

   if ((mrisimArgs::phantomFile != NULL) &&
       (mrisimArgs::fuzzy_specifier_string != NULL)){ 
      cerr << "Only one of -discrete or -fuzzy can be specified at a time."
           << endl;
   }
   if (mrisimArgs::tissueFile == NULL){
      cerr << "No tissue parameter file specified." << endl;
      exit(EXIT_FAILURE);
   }
   if (mrisimArgs::coilFile == NULL){
      cerr << "No coil parameter file specified." << endl;
      exit(EXIT_FAILURE);
   }
   if (mrisimArgs::sequenceFile == NULL){
      cerr << "No pulse sequence parameter file specified." << endl;
      exit(EXIT_FAILURE);
   }

   if (mrisimArgs::logFile != NULL){
      mrisimArgs::logFlag = TRUE;
   } else {
      mrisimArgs::logFlag = FALSE;
   }

   int count;
   if (mrisimArgs::fuzzy_specifier_string != NULL) {
      if (mrisimArgs::fuzzy_specifier_string[0] != '\0') {
         mrisimArgs::fuzzy_list = new Fuzzy_Specifier;
      }
      if (!parse_fuzzy_list(mrisimArgs::fuzzy_specifier_string, count, 
                            mrisimArgs::fuzzy_list)) {
         cerr << "Fuzzy phantom list parse failed after list element "
              << count << endl;
         delete_specifier_list(mrisimArgs::fuzzy_list); 
         mrisimArgs::fuzzy_list = NULL;
      }
   }

}

//------------------------------------------------------------------------- 
// mrisimArgs destructor
//------------------------------------------------------------------------- 

mrisimArgs::~mrisimArgs() {
   delete_fuzzy_list();
}

//------------------------------------------------------------------------- 
// mrisimArgs::get_phantom_type
// Returns the type of Phantom to create.
//------------------------------------------------------------------------- 

Phantom_Type mrisimArgs::get_phantom_type(void) const {
   Phantom_Type phantom_type;

   if (uses_discrete_phantom()) {
      if (uses_rf_phantom()) {
         phantom_type = DISCRETE_RF;
      } else {
         phantom_type = DISCRETE;
      }
   } else if (uses_fuzzy_phantom()) {
      if (uses_rf_phantom()) {
         phantom_type = FUZZY_RF;
      } else {
         phantom_type = FUZZY;
      }
   }
   return phantom_type;
}

//------------------------------------------------------------------------- 
// mrisimArgs::parse_optional_string_argument
// Handles parsing of an argv switch with an optional string argument.
// If the optional string argument is omitted, an empty string ('\0') is
// returned in dst, otherwise the string is copied to dst.
// If the string is copied, returns 1 to ParseArgv to mark the next
// argument for removal.  Otherwise, returns 0 to ParseArgv.
//------------------------------------------------------------------------- 

int parse_optional_string_argument(char *dst, char */*key*/, 
                                   char *nextArg) {

   int next_arg_used;

   // ParseArgv uses (char *) as a generic pointer.  This function is
   // always passed a (char **) as its first argument (cast to char *)
   // so some funky casting is required.

   if ((nextArg == NULL) || (nextArg[0] == '-')) {
      *(char **)dst = new char;
      (*(char **)dst)[0] = '\0';
      next_arg_used = 0;
   } else {
      *(char **)dst = strdup(nextArg);
      next_arg_used = 1;
   }

   return next_arg_used;
}

//------------------------------------------------------------------------- 
// mrisimArgs::delete_specifier_list
// Deletes a linked list of fuzzy phantom specifiers.
//------------------------------------------------------------------------- 

void mrisimArgs::delete_specifier_list(Fuzzy_Specifier *head) {
   Fuzzy_Specifier *ptr = head;

   while (head != NULL) {
      ptr = head->next;
      delete head;
      head = ptr;
   }
}

//------------------------------------------------------------------------- 
// mrisimArgs::split_string
// Splits a string into two parts, before and after the first occurence
// of the delimiter character.
// If the delimiter is found, returns TRUE, replaces the delimiter with
// '\0' and sets second_string to point to the character immediately
// following the delimiter.
// If the delimiter was not found, returns FALSE, and set second_string
// to point to the end-of-string character '\0'. 
//------------------------------------------------------------------------- 

int mrisimArgs::split_string(char *string, char delimiter, 
                             char **second_string) {
   int  found = FALSE;
   char ch;

   while ( ((ch=*string) != '\0') && !found ) {
      if (ch == delimiter)
         found = TRUE;
      else
         string++;
   }

   if (found) {
      *second_string = string+1;
      *string        = '\0';
   } else {
      *second_string = string;
   }

   return found;

}

//------------------------------------------------------------------------- 
// mrisimArgs::scan_label
// Validates an integer label by checking that each character in string
// is a digit. 
// If no non-digit characters were encountered, returns TRUE with label
// containing the converted integer label.
// If a non-digit character is encountered, returns FALSE with label
// containing the label formed from the valid digit characters.
//------------------------------------------------------------------------- 

int mrisimArgs::scan_label(const char *string, int &label) {
   const char *cptr = string;
 
   label = 0;
   while ( (*cptr != '\0') && (isdigit(*cptr)) ) {
      label = (*cptr-'0') + 10*label;
      cptr++;
   }

   return (*cptr == '\0' && cptr != string);
}

//------------------------------------------------------------------------- 
// mrisimArgs::parse_fuzzy_list
// Parse the fuzzy specifier list.  The list is contained in the string 
// head.  It consists of a zero-terminated string of comma-delimited
// tokens.  e.g.  <name-1>,<name-2>,...<name-N>.
// Optionally, each token also contains an integer label separated from
// the name by a colon.  e.g. <name-1>:<label-1>,...,<name-N>:<label-N>.
// The names and labels are parsed out and returned in the Fuzzy_Specifier
// linked list.   If labels are not used, the labels are set to -1.
// If a parse error occurs an error message is written and FALSE is
// returned, otherwise TRUE is returned.
//------------------------------------------------------------------------- 

int mrisimArgs::parse_fuzzy_list(char *head, 
                                 int &count, Fuzzy_Specifier *list) {
   int  more_members, uses_labels, has_label, no_error = TRUE;
   int  label;
   char *tail, *labelstr;

   count = 0;
   more_members = split_string(head, ',', &tail);
   uses_labels  = split_string(head, ':', &labelstr);

   if (*head != '\0') {
      if (uses_labels) {
         if (!scan_label(labelstr,label)) {
            cerr << "Improper label." << endl;
            no_error = FALSE;
         } else {
            //cout << "Name: " << head << " Label: " << label << endl;
            list->filename = head;
            list->label    = label;
            list->next     = NULL;
            count++;
         }
      } else {
         //cout << "Name: " << head << endl;
         list->filename = head;
         list->label    = -1;
         list->next     = NULL;
         count++;
      }
   }

   while (more_members && no_error) {
      head = tail;
      more_members = split_string(head, ',', &tail);
      if (*head == '\0') {
         cerr << "Empty list member not allowed." << endl;
         no_error = FALSE;
      } else {
         has_label = split_string(head, ':', &labelstr);
         if (uses_labels) {
            if (has_label) {
               if (!scan_label(labelstr,label)) {
                  cerr << "Improper label specified." << endl;
                  no_error = FALSE;
               } else {
                  //cout << "Name: " << head << " Label: " << label << endl;
                  list->next = new Fuzzy_Specifier;
                  list       = list->next;

                  list->filename = head;
                  list->label    = label;
                  list->next     = NULL;
                  count++;
               }
            } else {
               cerr << "Label specifier ':' expected." << endl;
               no_error = FALSE;
            }
         } else {
            if (has_label) {
               cerr << "Label specifier ':' not expected." << endl;
               no_error = FALSE;
            } else {
               //cout << "Name: " << head << endl;
               list->next = new Fuzzy_Specifier;
               list       = list->next;

               list->filename = head;
               list->label    = -1;
               list->next     = NULL;
               count++;
            }
         }
      }
   }

   return no_error;

}

