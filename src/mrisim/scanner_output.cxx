//==========================================================================
// SCANNER_OUTPUT.CXX
// Scanner_Output class.
// Inherits from:  
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995, 1996 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/scanner_output.cxx,v 1.4 2008-11-06 10:58:23 rotor Exp $
// $Log: scanner_output.cxx,v $
// Revision 1.4  2008-11-06 10:58:23  rotor
//  * fixed includes for iostream and friends
//  * updated for new release (1.0.2)
//
// Revision 1.3  2004/08/10 15:42:19  bert
// Fix a warning and avoid use of nonstandard ios::nocreate bit
//
// Revision 1.2  2003/06/11 11:36:05  crisco
//
// bug fix: the phase output image should not be multiplied by the gain!
//
// this was not a problem for Brainweb, because the 'phase' output type
// is not currently offered as an option ... but from now on it could be.
//
// Revision 1.1  2003/05/30 16:43:12  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  19:32:35  rkwan
// Release 3.1 - initial release.
//
// Revision 3.1  1996/07/19  19:31:13  rkwan
// Handles output file generation.
//
// Revision 3.1  1996/07/19  16:09:15  rkwan
// Release 3.1 - initial release.
//
//==========================================================================

#include <fstream>
#include "scanner_output.h"

//--------------------------------------------------------------------------
// Scanner_Output constructor
//--------------------------------------------------------------------------

Scanner_Output::Scanner_Output(const mrisimArgs &args,
                               char *time_stamp,
                               MRI_Scanner &scanner) 
   : _image(scanner.get_nrows(), scanner.get_ncols()) {

   // _good is set to FALSE if any file creation functions fail
   // in order to signal that the Scanner_Output object is not usable.
   _good = TRUE;

   // Set up the image type
   switch(scanner.get_image_type()){
      case REAL_IMAGE:
         _output_type = IMAGE_R;
         break;
      case IMAG_IMAGE:
         _output_type = IMAGE_I;
         break;
      case MOD_IMAGE:
         _output_type = IMAGE_M;
         break;
      case PHASE_IMAGE:
         _output_type = IMAGE_P;
         break;
   }

   // Set up reconstructed output files
   if (!create_output_file(args.outputFile, args.clobberFlag, 
                           time_stamp, scanner, _output)){
      _good = FALSE;
   }

   // Set up raw data files
   if (!args.oldpvFlag && scanner.save_raw_data()) {
      char *real_file_name = extend_path(args.outputFile,4,".raw_real");
      char *imag_file_name = extend_path(args.outputFile,4,".raw_imag");
      if (!create_output_file(real_file_name, args.clobberFlag,
                              time_stamp, scanner, _real_raw_data)){
         _good = FALSE;
      }
      if (!create_output_file(imag_file_name, args.clobberFlag,
                              time_stamp, scanner, _imag_raw_data)){
         _good = FALSE;
      }
      delete[] real_file_name;
      delete[] imag_file_name;
   }

}

//--------------------------------------------------------------------------
// Scanner_Output destructor
//--------------------------------------------------------------------------

Scanner_Output::~Scanner_Output() {}

//--------------------------------------------------------------------------
// Scanner_Output::create_output_file
// Creates and initializes an output MINC file.  Returns FALSE if
// file creation fails.
//--------------------------------------------------------------------------

int Scanner_Output::create_output_file(const char *path, int clobber,
                       const char *time_stamp, const MRI_Scanner &scanner,
                       O_MINC_File &minc_file) {

   // If clobber switch is not set, check to see if the file
   // already exists before creating it.

   if (!clobber && file_exists(path)){
      cerr << endl << "FATAL ERROR: Output file " << (char *)path 
           << " already exists." << endl;
      cerr << "Use a new file name or -clobber." << flush << endl;
      return FALSE;
   } 

   if (minc_file.create(path,(clobber ? NC_CLOBBER : NC_NOCLOBBER)) == 
          MI_ERROR) {
      cerr << endl << "FATAL ERROR: MINC could not create the output file " 
           << (char *)path << flush << endl;
      return FALSE;
   }

   // Set up the output file volume info and setup and attach an ICV

   if (!scanner.initialize_output(minc_file, time_stamp)) {
      cerr << endl << "FATAL ERROR: MINC output file " 
           << (char *)path << " setup failed." << flush << endl;
      return FALSE;
   }

   return TRUE;
}

//--------------------------------------------------------------------------
// Scanner_Output::display_info
// Displays info about the simulation and output volumes.
// Info is displayed to stdout or an optional log file.
//--------------------------------------------------------------------------

void Scanner_Output::display_info(const mrisimArgs &args,
                                  char *time_stamp, 
                                  const MRI_Scanner &scanner) const {

   if (args.verboseFlag)
      scanner.display_info(cout);
      cout << "Creating Output MINC file:" << endl;
      cout << "--------------------------" << endl << endl;
      _output.display_volume_info(cout);
   if (args.logFlag){
      ofstream logfile(((args.logFile[0] == '\0') ? 
                       "mrisim.log" : args.logFile));
      logfile << time_stamp << endl << endl;
      scanner.display_info(logfile);
      logfile << "Creating Output MINC file:" << endl;
      logfile << "--------------------------" << endl << endl;
      _output.display_volume_info(logfile);
      logfile.close();
   }

}

//--------------------------------------------------------------------------
// Scanner_Output::save_images
// Generates simulated images from pulse sequence simulation and phantom
// data, and writes the simulated images in output files.
//--------------------------------------------------------------------------

void Scanner_Output::save_images(const mrisimArgs &args, 
                                 MRI_Scanner &scanner) {

   if (args.verboseFlag)
      cout << "Saving slices";

   int islice;

   if (args.oldpvFlag) {

      for(islice=0; islice<_output.get_nslices(); islice++){
         scanner.get_simulated_image_slice(islice, _image);
         _image.scale(scanner.get_signal_gain());
         _output.save_slice(islice, _image);
         if (args.verboseFlag)
            cout << "." << flush;
      }

   } else {

      scanner.initialize_chirp_resample();
      Complex_Slice raw_slice(scanner.get_matrix_size(ROW), 
                              scanner.get_matrix_size(COLUMN));

      for(islice=0; islice<_output.get_nslices(); islice++){
         scanner.get_raw_data_slice(islice, raw_slice);

         if (scanner.save_raw_data()){
            scanner.get_real_image(raw_slice, _image);
            _real_raw_data.save_slice(islice, _image);
            scanner.get_imag_image(raw_slice, _image);
            _imag_raw_data.save_slice(islice, _image);
         }

         // --- Save reconstructed image --- //

         scanner.reconstruct_raw_data_slice(raw_slice);
         switch(_output_type) {
            case IMAGE_R:
               scanner.get_real_image(raw_slice, _image);
               break;
            case IMAGE_I:
               scanner.get_imag_image(raw_slice, _image);
               break;
            case IMAGE_P:
               scanner.get_angle_image(raw_slice, _image);
               break;
            case IMAGE_M:
            default:
               scanner.get_abs_image(raw_slice, _image);
               break;
         }

	 // [CC, 11.06.2003] 
	 // multiplication with gain doesn't make sense for the phase 
	 if (_output_type != IMAGE_P) 
	   _image.scale(scanner.get_signal_gain());

         _output.save_slice(islice, _image);
         if (args.verboseFlag)
            cout << "." << flush;
      }

   }

   cout << endl;

}

//--------------------------------------------------------------------------
// Scanner_Output::file_exists
// Returns TRUE if a file exists.
//--------------------------------------------------------------------------

int Scanner_Output::file_exists(const char *path) const {
   int exist;
   ofstream tmp(path, ios::in);
   if (!tmp){
      exist = 0;
   } else {
      exist = 1;
   }
   tmp.close();
   return exist;

}

//--------------------------------------------------------------------------
// Scanner_Output::extend_path
// Returns a new string, with extend_string inserted in original_path
// extension_length characters from the end (not included terminating
// zero).
//--------------------------------------------------------------------------

char *Scanner_Output::extend_path(const char *original_path, 
                  int extension_length,
                  const char *extend_string) const {

   const int pathlen = strlen(original_path);
   const int extlen  = strlen(extend_string);

   char *extended_path = new char[pathlen+extlen+1]; // include space for '\0'

   memcpy(extended_path, 
          original_path, 
          (pathlen-extension_length)*sizeof(char));
   memcpy(&extended_path[pathlen-extension_length],
          extend_string, 
          extlen*sizeof(char));
   memcpy(&extended_path[pathlen-extension_length+extlen],
          &original_path[pathlen-extension_length],
          extension_length*sizeof(char));
   extended_path[pathlen+extlen] = '\0';

   return extended_path;

}

