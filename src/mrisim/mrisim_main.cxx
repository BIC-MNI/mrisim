//==========================================================================
// MRISIM_MAIN.CXX
// Main file for mrisim command line version.
// Inherits from:
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/mrisim_main.cxx,v 1.3 2008-11-06 10:58:23 rotor Exp $
// $Log: mrisim_main.cxx,v $
// Revision 1.4  2017-04-03 16:00:00 llewis
//  * added hires parameter allowances
//
// Revision 1.3  2008-11-06 10:58:23  rotor
//  * fixed includes for iostream and friends
//  * updated for new release (1.0.2)
//
// Revision 1.2  2004/08/10 15:37:38  bert
// Fix two bugs and a warning
//
// Revision 1.1  2003/05/30 16:43:11  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  15:58:33  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:15:42  rkwan
// Release 2.5
//
// Revision 1.5  1996/01/17  18:02:01  rkwan
// Update for tx/rx coil modelling.
//
// Revision 1.4  1995/12/22  20:16:57  rkwan
// Update for percent coil and RF map features.
//
// Revision 1.3  1995/12/11  17:50:40  rkwan
// Check for output file existence added.
//
// Revision 1.2  1995/12/11  15:24:43  rkwan
// Doc update.
//
//==========================================================================

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string.h>
#include <ctype.h>

#include "mrisim_main.h"
#include "scanner_output.h"

//--------------------------------------------------------------------------
// main
//--------------------------------------------------------------------------

int main(int argc, char *argv[]){

   // --- INITIALIZATION --- //

   // Save the time stamp
   char *stamp = time_stamp(argc, argv);

   // Check command line arguments
   mrisimArgs args(argc, argv);

   // --- CREATE SIMULATOR MODELS --- //

   MRI_Scanner  scanner;
   create_models(args, scanner);
   //args.delete_fuzzy_list();

   // --- APPLY PULSE SEQUENCE --- //
   // Create a Pulse_Sequence from the parameter file and apply
   // it to the MRI_Scanner model.

   if (!apply_pulse_sequence(args, scanner)){
       cerr << endl << "FATAL ERROR: Could not create Pulse Sequence."
            << flush << endl;
       exit(EXIT_FAILURE);
   }

   // --- OUTPUT IMAGES --- //

   Scanner_Output output(args, stamp, scanner);

   if (!output.is_good()) {
      exit(EXIT_FAILURE);
   }

   output.display_info(args, stamp, scanner);
   output.save_images(args, scanner);

   // --- CLEAN UP --- //
   free(stamp);

   return 0;

}

//--------------------------------------------------------------------------
// create_models
//
// Initialize the MRI_Scanner model.
// Parse parameter files and create Phantom and RF_Coil models, attach
// these models to the MRI_Scanner.   These models become the property
// of the MRI_Scanner and are destroyed automatically with the 
// MRI_Scanner.
//--------------------------------------------------------------------------

void create_models(const mrisimArgs &args, MRI_Scanner &scanner) {

   // --- SCANNER INITIALIZATION --- //

   // Set the scanner signal gain
   scanner.set_signal_gain(args.signal_gain);

   // --- RF_COIL INITIALIZATION --- //
   
   // Create an RF_Coil model from parameter files and attach it
   // to the MRI_Scanner

   RF_Coil *rf_coil;
   if ((rf_coil = make_coil(args)) == NULL){
       cerr << endl << "FATAL ERROR: Could not create RF coil."
            << flush << endl;
       exit(EXIT_FAILURE);
   }
   if (!scanner.attach(rf_coil)){
       cerr << endl << "FATAL ERROR: Could not attach RF coil."
            << flush << endl;
       exit(EXIT_FAILURE);
   }
   
   // --- PHANTOM INITIALIZATION --- //

   // Create a Phantom model from parameter files and attach it
   // to the MRI_Scanner

   Phantom *phantom;
   if ((phantom = open_phantom(args, rf_coil)) == NULL){
      cerr << endl << "FATAL ERROR: Could not create phantom." 
           << flush << endl;
      exit(EXIT_FAILURE);
   }
   if (!scanner.attach(phantom)){
      cerr << endl << "FATAL ERROR: Could not attach phantom."
           << flush << endl;
      exit(EXIT_FAILURE);
   }

}

//--------------------------------------------------------------------------
// open_phantom
// Parses the mrisim phantom and tissue parameter files and creates a new
// Phantom model.
// Returns FALSE if phantom creation has failed.
//--------------------------------------------------------------------------

Phantom *open_phantom(const mrisimArgs& args, const RF_Coil *rf_coil) {

   int       status = TRUE;

   // --- Phantom Pointers --- //

   Discrete_Phantom    *dphantom   = NULL;
   Discrete_RF_Phantom *drfphantom = NULL;
   Fuzzy_Phantom       *fphantom   = NULL;
   Fuzzy_RF_Phantom    *frfphantom = NULL;
   Phantom             *phantom    = NULL;

   ParamFile paramfile(args.tissueFile);

   // --- Input Parameters --- //

   int  n_tissue_classes, n_flip_angles;
   char phantom_type_specifier[255]; 
   int  max_tissue_label;
   char label_file_name[255];

   char         tissue_name[255];
   int          tissue_label;
   double       T1, T2, T2s, NH;

   // --- Read phantom tissue description parameters --- //

   paramfile.getfield(n_tissue_classes);
   if (n_tissue_classes < 1 || n_tissue_classes > 255) {
      cerr << "Number of tissue classes " << n_tissue_classes 
           << " is out of range 1 - 255." << endl;
      status = FALSE;
   }

   paramfile.getfield(n_flip_angles);
   if (n_flip_angles < 1 || n_flip_angles > 10000) {
      cerr << "Number of flip angles " << n_flip_angles
           << " is out of range 1 - 10000." << endl;
      status = FALSE;
   }

   // If -discrete or -fuzzy were not specified on the command-line
   // use values from the parameter file.

   paramfile.getfield(255, phantom_type_specifier);
   if (!(args.uses_fuzzy_phantom() || args.uses_discrete_phantom())) {
      if (strcmp(phantom_type_specifier, "discrete")==0) {
         args.phantomFile = new char;
         args.phantomFile[0] = '\0';
      } else if (strcmp(phantom_type_specifier, "fuzzy")==0) {
         args.fuzzy_specifier_string = new char;
         args.fuzzy_specifier_string[0] = '\0';
      } else {
         cerr << "Invalid phantom type: " << phantom_type_specifier << endl;
         cerr << "Phantom Type must be one of: discrete, fuzzy."
              << endl;
         status = FALSE;
      }
   }

   // Determine the phantom type
   Phantom_Type phantom_type = args.get_phantom_type();

   paramfile.getfield(max_tissue_label);
   if (max_tissue_label < 0 || max_tissue_label > 255) {
      cerr << "Highest tissue label " << max_tissue_label
           << " is out of range 0 - 255." << endl;   
      status = FALSE;
   }

   // Read the default discrete phantom file name.
   // If args.uses_default_discrete_phantom() is TRUE this
   // file name will be used.
   paramfile.getfield(255, label_file_name);

   // --- Check for file parsing errors --- //
 
   if (status == FALSE) {
      return NULL;
   }

   // --- Create Phantom Instances --- //

   switch (phantom_type) {
      case DISCRETE:      // --- Discrete Phantom --- //

         dphantom   = new Discrete_Phantom(n_tissue_classes);
         if (args.uses_default_discrete_phantom()) {
            dphantom->open_discrete_label_file(label_file_name);
         } else {
            dphantom->open_discrete_label_file(args.phantomFile);
         }
         phantom    = (Phantom *)dphantom;
         break;

      case DISCRETE_RF:   // --- Discrete RF Phantom --- //

         drfphantom = new Discrete_RF_Phantom(n_tissue_classes, n_flip_angles);
         if (args.uses_default_discrete_phantom()) {
            drfphantom->open_discrete_label_file(label_file_name);
         } else {
            drfphantom->open_discrete_label_file(args.phantomFile);
         }
         phantom    = (Phantom *)drfphantom;
         break;

      case FUZZY:       // --- Fuzzy Phantom --- //

         fphantom   = new Fuzzy_Phantom(n_tissue_classes); 
         phantom    = (Phantom *)fphantom;
         break;

      case FUZZY_RF:    // --- Fuzzy RF Phantom --- //
   
         frfphantom = new Fuzzy_RF_Phantom(n_tissue_classes, n_flip_angles);
         phantom    = (Phantom *)frfphantom;
         break;
   }

   // --- Read tissue parameters and install tissues in the phantom --- //

   Fuzzy_Specifier *fuzzy_list = args.fuzzy_list;

   while (n_tissue_classes-- > 0) {

      paramfile.getfield(255, tissue_name);

      paramfile.getfield(tissue_label);
      if (tissue_label > 255){
         cerr << "Tissue Label " << tissue_label 
              << " is out of range 0 - 255" << endl;
         status = FALSE;
      }

      paramfile.getfield(255, label_file_name);

      paramfile.getfield(T1);
      if (T1 < 0.0 || T1 > 20000.0){
         cerr << "T1 " << T1 << " for tissue " << tissue_label
              << " is out of range 0.00 - 20000.00" << endl;
         status = FALSE;
      }
 
      paramfile.getfield(T2);
      if (T2 < 0.0 || T2 > 20000.0){
         cerr << "T2 " << T2 << " for tissue " << tissue_label
              << " is out of range 0.00 - 20000.00" << endl;
         status = FALSE;
      }

      paramfile.getfield(T2s);
      if (T2s < 0.0 || T2s > 20000.0){
         cerr << "T2* " << T2s << " for tissue " << tissue_label
              << " is out of range 0.00 - 20000.00" << endl;
         status = FALSE;
      }

      paramfile.getfield(NH);
      if (NH < 0.0 || NH > 1.0) {
         cerr << "Proton Density " << NH << " for tissue " << tissue_label
              << " is out of range 0.00 - 1.00" << endl;
         status = FALSE;
      }

      // Do not install tissues with 0 proton density if they are not
      // background (tissue_label of 0)
      if (NH != 0.0 || tissue_label == 0 || tissue_label == 7){

         // Install tissue parameters in phantom
         phantom->install_tissue(new Tissue(tissue_name,T1,T2,T2s,NH),
                                 tissue_label);

         // If fuzzy labels were not specified store them now
         //if (args.uses_fuzzy_phantom() && !args.fuzzy_labels_specified()) {
         //   fuzzy_list->label = tissue_label;
         //   fuzzy_list = fuzzy_list->next;
         //}
        
         // If a fuzzy phantom is used with default label files
         // then open the files now.
         if (args.uses_default_fuzzy_phantom()) {

            switch (phantom_type) {
               case FUZZY:
                  if (!fphantom->open_fuzzy_label_file(tissue_label,
                                    label_file_name)) {
                     cerr << "Invalid fuzzy label file: " 
                          << label_file_name << " for tissue label " 
                          << tissue_label << "." << endl;
                     status = FALSE;
                  }
                  break;
               case FUZZY_RF:
                  if (!frfphantom->open_fuzzy_label_file(tissue_label,
                                     label_file_name)) {
                     cerr << "Invalid fuzzy label file: "
                          << label_file_name << " for tissue label " 
                          << tissue_label << "." << endl;
                     status = FALSE;
                  }
                  break;
               default:
                  break;
            }
         }
      }
 
      if (status == FALSE) {
         delete phantom;
         return NULL;
      }

   }

   if (args.uses_fuzzy_phantom() && !args.uses_default_fuzzy_phantom()) {
      fuzzy_list = args.fuzzy_list;
      while(fuzzy_list != NULL) {
         // Check if the tissue label is valid.
         if (!phantom->tissue_label_is_valid(fuzzy_list->label)) {
            cerr << "Error in fuzzy phantom specification." << endl;
            cerr << "No tissue label " << fuzzy_list->label 
                 << " specified in the tissue parameter file."
                 << endl;
            delete phantom;
            return NULL;
         }
         switch (phantom_type) {
            case FUZZY:
               if (!fphantom->open_fuzzy_label_file(fuzzy_list->label,
                                 fuzzy_list->filename)) {
                  cerr << "Invalid fuzzy label file: " 
                       << fuzzy_list->filename << " for tissue label " 
                       << fuzzy_list->label << "." << endl;
                  status = FALSE;
               }
               break;
            case FUZZY_RF:
               if (!frfphantom->open_fuzzy_label_file(fuzzy_list->label,
                                  fuzzy_list->filename)) {
                  cerr << "Invalid fuzzy label file: "
                       << fuzzy_list->filename << " for tissue label " 
                       << fuzzy_list->label << "." << endl;
                  status = FALSE;
               }
               break;
            default:
               break;
         }
      
         if (status == FALSE) {
            delete phantom;
            return NULL;
         }

         fuzzy_list = fuzzy_list->next;
      }
   }

   // --- RF Map Consistency Checks --- //

   switch (phantom_type) {
      case DISCRETE_RF:
         if (!drfphantom->use_rf_coil(rf_coil)) {
            delete drfphantom;
            return NULL;
         }
         if (rf_coil->uses_rx_map() && 
             !drfphantom->has_same_dimension_info_as_rx(*rf_coil)){
            cerr << "WARNING: Receive map dimension info not the "
                 << "same as Phantom." << flush << endl;
         }
         if (rf_coil->uses_tx_map() && 
             !drfphantom->has_same_dimension_info_as_tx(*rf_coil)){
            cerr << "WARNING: Transmit map dimension info not the "
                 << "same as Phantom." << flush << endl;
            }
         if (rf_coil->uses_rx_map() && 
             !drfphantom->has_same_dimension_names_as_rx(*rf_coil)){
            cerr << "WARNING: Receive map dimension names not the "
                 << "same as Phantom." << flush << endl;
         }
         if (rf_coil->uses_tx_map() && 
             !drfphantom->has_same_dimension_info_as_tx(*rf_coil)){
            cerr << "WARNING: Transmit map dimension info not the "
                 << "same as Phantom." << flush << endl;
         }
         break;
      case FUZZY_RF:
         if (!frfphantom->use_rf_coil(rf_coil)) {
            delete frfphantom;
            return NULL;
         }
         if (rf_coil->uses_rx_map() && 
             !frfphantom->has_same_dimension_info_as_rx(*rf_coil)){
            cerr << "WARNING: Receive map dimension info not the "
                 << "same as Phantom." << flush << endl;
         }
         if (rf_coil->uses_tx_map() && 
             !frfphantom->has_same_dimension_info_as_tx(*rf_coil)){
            cerr << "WARNING: Transmit map dimension info not the "
                 << "same as Phantom." << flush << endl;
            }
         if (rf_coil->uses_rx_map() && 
             !frfphantom->has_same_dimension_names_as_rx(*rf_coil)){
            cerr << "WARNING: Receive map dimension names not the "
                 << "same as Phantom." << flush << endl;
         }
         if (rf_coil->uses_tx_map() && 
             !frfphantom->has_same_dimension_info_as_tx(*rf_coil)){
            cerr << "WARNING: Transmit map dimension info not the "
                 << "same as Phantom." << flush << endl;
         }
         break;
   }

   paramfile.close();

   return phantom;

}

//--------------------------------------------------------------------------
// make_coil
// Parses the mrisim rf_coil parameter file and creates a new coil
// model.
//--------------------------------------------------------------------------

RF_Coil *make_coil(const mrisimArgs& args) {

   RF_Coil   *rf_coil = NULL;

   int       status = TRUE;
   ParamFile paramfile(args.coilFile);

   enum Noise_Model {NOISELESS=0, PERCENT=1, INTRINSIC=2, IMAGESNR=3};
   Noise_Model noise_model;
 
   double   intrinsic_snr, percent_noise, image_snr;
   long     random_seed;
   char     noise_model_specifier[255];
   double   reference_thickness = 0;
   int      reference_tissue = 0;
   char     rxmap[255];
   char     txmap[255];
 
   // --- Read coil parameter description --- //

   paramfile.getfield(255, noise_model_specifier);
   if (strcmp(noise_model_specifier, "noiseless")==0){
      noise_model = NOISELESS;
   } else if (strcmp(noise_model_specifier, "percent")==0){
      noise_model = PERCENT;
   } else if (strcmp(noise_model_specifier, "intrinsic")==0){
      noise_model = INTRINSIC;
   } else if (strcmp(noise_model_specifier, "imagesnr")==0){
      noise_model = IMAGESNR;
   } else {
      cerr << "Invalid noise model: " << noise_model_specifier << "." << endl;
      cerr << "Noise model must be one of: "  
           << "noiseless, percent, intrinsic, imagesnr" << endl;
      status = FALSE;
   }
 
   paramfile.getfield(percent_noise);
   if (percent_noise < 0.0 || percent_noise > 100.0){
      cerr << "Percent noise " << percent_noise 
           << " is out of range 0.00 - 100.00" << endl; 
      status = FALSE;
   }

   paramfile.getfield(reference_thickness);
   if ((reference_thickness < 0.0001 || reference_thickness > 320.00) &&
       (reference_thickness != 0.0))  {
      cerr << "Reference slice thickness " << reference_thickness
           << " is out of range 0.0001 - 320.00" << endl;
      status = FALSE;
   }

   paramfile.getfield(reference_tissue);
   if ((reference_tissue < 1 || reference_tissue > 255) &&
       (reference_tissue != 0)) {
      cerr << "Reference tissue " << reference_tissue
           << " is out of range 1 - 255" << endl;
      status = FALSE;
   }

   paramfile.getfield(intrinsic_snr);
   paramfile.getfield(255, rxmap);
   if (args.uses_rx_map()) {
      if (args.uses_default_rx_map()) { // args.rxmapFile[0]=='\0'
         if (rxmap[0] != '\0') {
            delete[] args.rxmapFile;
            args.rxmapFile = strdup(rxmap);
         } else {
            cerr << "No default receive map specified. "
                 << "Ignoring -rxmap switch." << endl;
            delete[] args.rxmapFile;
            args.rxmapFile = NULL;
         }
      }
   } else if (rxmap[0] != '\0') {
      args.rxmapFile = strdup(rxmap);
   }

   paramfile.getfield(255, txmap);
   if (args.uses_tx_map()) {
      if (args.uses_default_tx_map()) { // args.txmapFile[0]=='\0'
         if (txmap[0] != '\0') {
            delete[] args.txmapFile;
            args.txmapFile = strdup(txmap);
         } else {
            cerr << "No default transmit map specified. "
                 << "Ignoring -txmap switch." << endl;
            delete[] args.txmapFile;
            args.txmapFile = NULL;
         }
      }
   } else if (txmap[0] != '\0') {
      args.txmapFile = strdup(txmap);
   }

   paramfile.getfield(random_seed);

   // -- Check for file parsing errors --- //

   if (status == FALSE){
      return NULL;
   }
  
   // --- Generate RF Coil Instances --- //
 
   switch(noise_model){
      case NOISELESS:
         rf_coil = new RF_Coil(random_seed, args.rxmapFile, args.txmapFile);
         break;
      case PERCENT:
         rf_coil = new Percent_Coil(percent_noise, reference_thickness, 
                                    reference_tissue, random_seed, 
                                    args.rxmapFile, args.txmapFile);
         break;
      case IMAGESNR:
         rf_coil = new Image_SNR_Coil(image_snr, random_seed, 
                                    args.rxmapFile, args.txmapFile);
         break;
      default:
      case INTRINSIC:
         rf_coil = new Intrinsic_SNR_Coil(intrinsic_snr, random_seed, 
                                    args.rxmapFile, args.txmapFile);
         break;
   }
   paramfile.close();

   return rf_coil;
}

//--------------------------------------------------------------------------
// apply_pulse_sequence
// Parses the mrisim pulse sequence parameter file and creates a
// new pulse sequence model.
//--------------------------------------------------------------------------

int apply_pulse_sequence(const mrisimArgs &args, MRI_Scanner &scanner) {

   int       status = TRUE;
   ParamFile paramfile(args.sequenceFile);

   // --- Pulse Sequence Pointers --- //
   Quick_Sequence  *quick_pseq;
   Custom_Sequence *custom_pseq;
   Pulse_Sequence  *pseq;

   // --- Pulse Sequence Parameters --- //
   double    voxel_offset[3];
   int       foldover_direction, foldover_suppression;
   int       number_of_slices;
   double    slice_thickness, slice_separation, field_of_view;
   double    rectangular_fov, scan_percentage;
   double    TR, TI, TE1, TE2;
   int       number_of_echoes, partial_echo;
   double    flip_angle, water_fat_shift;
   int       number_of_signal_averages;
   int       half_scan, scan_matrix, reconstruction_matrix;
   char      buffer[255];
 
   enum Scan_Technique {SCAN_TYPE_SE=0, SCAN_TYPE_IR=1, SCAN_TYPE_SFLASH=2, 
                        SCAN_TYPE_CEFAST=3, SCAN_TYPE_FISP=4, SCAN_TYPE_FLASH=5,
                        SCAN_TYPE_DSE_EARLY=6, SCAN_TYPE_DSE_LATE=7};

   Scan_Technique    scan_technique;
   Image_Orientation orientation;
   Scan_Mode         scan_mode;
   Image_Type        image_type;
 
   Phantom   *phantom = scanner.get_attached_phantom();

   // --- Read in Pulse Sequence Parameter File --- //

   paramfile.getfield(voxel_offset[SLICE]);
   paramfile.getfield(voxel_offset[ROW]);
   paramfile.getfield(voxel_offset[COLUMN]);
   if (args.voxel_offset[SLICE] != DBL_MIN) {
      voxel_offset[SLICE] = args.voxel_offset[SLICE];
   }      
   if (args.voxel_offset[ROW] != DBL_MIN) {
      voxel_offset[ROW] = args.voxel_offset[ROW];
   }      
   if (args.voxel_offset[COLUMN] != DBL_MIN) {
      voxel_offset[COLUMN] = args.voxel_offset[COLUMN];
   }      

   if (voxel_offset[SLICE] < -160.0 || voxel_offset[SLICE] > 160.0){
      cerr << "Slice voxel offset " << voxel_offset[SLICE]
           << " is out of range: -160.00 - 160.00." << endl;
      status = FALSE;
   }
   if (voxel_offset[ROW] < -220.0 || voxel_offset[ROW] > 220.0){
      cerr << "Row voxel offset " << voxel_offset[ROW]
           << " is out of range: -220.00 - 220.00." << endl;
      status = FALSE;
   }
   if (voxel_offset[COLUMN] < -220.0 || voxel_offset[COLUMN] > 220.0){
      cerr << "Column voxel offset " << voxel_offset[COLUMN]
           << " is out of range: -220.00 - 220.00." << endl;
      status = FALSE;
   }
  
   paramfile.getfield(255, buffer); // ignore slice orientation
   orientation = SAME;
   paramfile.getfield(255, buffer);
   if (strcmp(buffer,"X")==0){
      foldover_direction = COLUMN;
   } else if (strcmp(buffer,"Y")==0){
      foldover_direction = ROW;
   } else {
      cerr << "Invalid foldover direction: " << buffer << endl;
      cerr << "Foldover direction must be one of: Y, X." << endl;
      status = FALSE;
   }

   paramfile.getfield(255, buffer);
   if (strcmp(buffer,"no")==0){
      foldover_suppression = NO;
   } else if (strcmp(buffer,"yes")==0){
      foldover_suppression = YES;
   } else {
      cerr << "Invalid foldover suppression: " << buffer << endl;
      cerr << "Foldover suppression must be one of: no, yes." << endl;
      status = FALSE;
   }

   paramfile.getfield(number_of_slices);
   if (number_of_slices < 0 || number_of_slices > 1000000) {
      cerr << "Number of slices " << number_of_slices
           << " is out of range: 1 - 1000000." << endl;
      status = FALSE;
   }

   paramfile.getfield(slice_thickness);
   if (slice_thickness < 0.0001 || slice_thickness > 320.0) {
      cerr << "Slice thickness " << slice_thickness
           << " is out of range: 0.0001 - 320.00." << endl;
      status = FALSE;
   }

   paramfile.getfield(slice_separation);
   if (slice_separation < 0.0001 || slice_separation > 320.0) {
      cerr << "Slice separation " << slice_separation
           << " is out of range: 0.0001 - 320.00." << endl;
      status = FALSE;
   }

   paramfile.getfield(field_of_view);
   if (field_of_view < 40.00 || field_of_view > 500.00){
      cerr << "Field of view " << field_of_view
           << " is out of range: 40.00 - 500.00." << endl;
      status = FALSE;
   }

   paramfile.getfield(rectangular_fov);
   if (rectangular_fov < 25.0 || rectangular_fov > 100.0){
      cerr << "Rectangular FOV " << rectangular_fov
           << " is out of range: 25.00 - 100.00." << endl;
      status = FALSE;
   }

   paramfile.getfield(255, buffer);
   if (strcmp(buffer,"SE")==0){
      scan_technique = SCAN_TYPE_SE;
   } else if (strcmp(buffer,"IR")==0){
      scan_technique = SCAN_TYPE_IR;
   } else if (strcmp(buffer,"SFLASH")==0){
      scan_technique = SCAN_TYPE_SFLASH;
   } else if (strcmp(buffer,"CEFAST")==0){
      scan_technique = SCAN_TYPE_CEFAST;
   } else if (strcmp(buffer,"FISP")==0){
      scan_technique = SCAN_TYPE_FISP;
   } else if (strcmp(buffer,"FLASH")==0){
      scan_technique = SCAN_TYPE_FLASH;
   } else if (strcmp(buffer,"DSE_EARLY")==0){
      scan_technique = SCAN_TYPE_DSE_EARLY;
   } else if (strcmp(buffer,"DSE_LATE")==0){
      scan_technique = SCAN_TYPE_DSE_LATE;
   } else {
      cerr << "Invalid scan technique: " << buffer << endl;
      cerr << "Scan technique must be one of: SE, IR, SFLASH, CEFAST, "
           << "FISP, FLASH, DSE_EARLY, DSE_LATE." << endl;
      status = FALSE;
   }

   paramfile.getfield(255, buffer);
   if (strcmp(buffer,"2D")==0){
      scan_mode = SCAN_MODE_2D;
   } else if (strcmp(buffer,"3D")==0){
      scan_mode = SCAN_MODE_3D;
   } else if (strcmp(buffer,"MS")==0){
      scan_mode = SCAN_MODE_MS;
   } else {
      cerr << "Invalid scan mode: " << buffer << endl;
      cerr << "Scan mode must be one of: 2D, 3D, MS." << endl;
      status = FALSE;
   }

   paramfile.getfield(TR);
   if (TR < 0.0) {
      cerr << "Repetition time " << TR << " must be positive." << endl;
      status = FALSE;
   }

   paramfile.getfield(TI);
   if (TI > TR) {
      cerr << "Inversion time " << TI 
           << " must be less than repetition time" << TR << "." << endl;
      status = FALSE;
   }

   paramfile.getfield(number_of_echoes);
   if (number_of_echoes < 1 || number_of_echoes > 10){
      cerr << "Number of echoes " << number_of_echoes
           << " is out of range: 1 - 10." << endl;
      status = FALSE;
   }

   paramfile.getfield(255, buffer);
   if (strcmp(buffer,"no")==0){
      partial_echo = NO;
   } else if (strcmp(buffer,"yes")==0){
      partial_echo = YES;
   } else {
      cerr << "Invalid partial echo: " << buffer << endl;
      cerr << "Partial echo must be one of: no, yes." << endl;
      status = FALSE;
   }

   paramfile.getlist2(TE1, TE2);
   paramfile.getfield(flip_angle);
   if (flip_angle < 1.0 || flip_angle > 150.0){
      cerr << "Flip angle " << flip_angle
           << " is out of range: 1.00 - 150.00." << endl;
      status = FALSE;
   }

   paramfile.getfield(water_fat_shift);
   if (water_fat_shift < 0.5 || water_fat_shift > 6.0){
      cerr << "Water fat shift " << water_fat_shift 
           << " is out of range: 0.50 - 6.00." << endl;
      status = FALSE;
   }

   paramfile.getfield(number_of_signal_averages);
   if (number_of_signal_averages < 1 || number_of_signal_averages > 32){
      cerr << "Number of signals averaged " << number_of_signal_averages
           << " is out of range: 1 - 32." << endl;
      status = FALSE;
   }

   paramfile.getfield(255, buffer);
   if (strcmp(buffer,"no")==0){
      half_scan = NO;
   } else if (strcmp(buffer,"yes")==0){
      half_scan = YES;
   } else {
      cerr << "Invalid half scan: " << buffer << endl;
      cerr << "Half scan must be one of: no, yes." << endl;
      status = FALSE;
   }

   paramfile.getfield(scan_percentage);
   if (scan_percentage < 25.0 || scan_percentage > 100.0){
      cerr << "Scan percentage " << scan_percentage
           << " is out of range: 25.00 - 100.00." << endl;
      status = FALSE;
   }

   paramfile.getfield(scan_matrix);
   if ((scan_matrix != 64)  && (scan_matrix != 128) && (scan_matrix != 256) && 
       (scan_matrix != 512) && (scan_matrix != 1024) && (scan_matrix != 2048) && (scan_matrix != 4096) && (scan_matrix != 8192)){   
      cerr << "Invalid scan matrix: " << scan_matrix << endl;
      cerr << "Scan matrix must be one of: 64, 128, 256, 512, 1024, 2048, 4096, 8192." << endl;
      status = FALSE;
   }

   paramfile.getfield(reconstruction_matrix);
   if ((reconstruction_matrix != 64)  && (reconstruction_matrix != 128) && (reconstruction_matrix != 256) &&
       (reconstruction_matrix != 512) && (reconstruction_matrix != 1024) && (reconstruction_matrix != 2048) && (reconstruction_matrix != 4096) && (reconstruction_matrix != 8192)){   
      cerr << "Invalid reconstruction matrix: " << reconstruction_matrix 
           << endl;
      cerr << "Reconstruction matrix must be one of: 64, 128, 256, 512, 1024, 2048, 4096, 8192." 
           << endl;
      status = FALSE;
   }

   paramfile.getfield(255, buffer);
   if (strcmp(buffer,"M")==0){
      image_type = MOD_IMAGE;
   } else if (strcmp(buffer,"R")==0){
      image_type = REAL_IMAGE;
   } else if (strcmp(buffer,"I")==0){
      image_type = IMAG_IMAGE;
   } else if (strcmp(buffer,"P")==0){
      image_type = PHASE_IMAGE;
   } else if (strcmp(buffer,"none")==0){
      image_type = NO_IMAGE;
   } else {
      cerr << "Invalid image type: " << buffer << endl;
      cerr << "Image type must be one of: R, I, M, P, none." << endl;
      status = FALSE;
   }

   paramfile.getfield(255, buffer);
   if (strcmp(buffer,"no")==0){
      scanner.save_raw_data(NO);
   } else if (strcmp(buffer,"yes")==0){
      scanner.save_raw_data(YES);
   } else {
      cerr << "Invalid save raw data parameter: " << buffer << endl;
      cerr << "Save raw data must be one of: no, yes." << endl;
      status = FALSE;
   }

   // --- Check for parameter input errors --- //
   if (status == FALSE) {
      return FALSE;
   }

   // --- Generate Pulse Sequence --- //

   RF_Pulse pulse1(X_AXIS,(Degrees)flip_angle);
   RF_Pulse pulse2(X_AXIS,(Degrees)180.0);
   RF_Pulse pulse3(X_AXIS,(Degrees)180.0);

   switch (scan_technique) {
      case SCAN_TYPE_SE:
         quick_pseq = new SE(TR,TE1,&Phantom::_save_sample, (void *)phantom);
         pseq = (Pulse_Sequence *)quick_pseq;
         break;
      case SCAN_TYPE_SFLASH:
         quick_pseq = 
            new Spoiled_FLASH(TR,TE1,flip_angle,&Phantom::_save_sample,
                                  (void *)phantom);
         pseq = (Pulse_Sequence *)quick_pseq;
         break;
      case SCAN_TYPE_IR:
         quick_pseq = 
            new IR(TR,TE1,TI,&Phantom::_save_sample,
                       (void *)phantom);
         pseq = (Pulse_Sequence *)quick_pseq;
         break;
      case SCAN_TYPE_CEFAST:
         quick_pseq = 
            new CE_FAST(TR,TE1,flip_angle,&Phantom::_save_sample,
                            (void *)phantom);
         pseq = (Pulse_Sequence *)quick_pseq;
         break;
      case SCAN_TYPE_FISP:
         quick_pseq = 
            new FISP(TR,TE1,flip_angle,&Phantom::_save_sample,
                         (void *)phantom);
         pseq = (Pulse_Sequence *)quick_pseq;
         break;
      case SCAN_TYPE_FLASH:
         quick_pseq = 
            new FLASH(TR,TE1,flip_angle,&Phantom::_save_sample,
                          (void *)phantom);
         pseq = (Pulse_Sequence *)quick_pseq;
         break;
      case SCAN_TYPE_DSE_EARLY:
         custom_pseq = new Custom_Sequence;
         if (TE2 < TE1){
            cerr << "First echo time must be less than second echo time." 
                 << endl;
            delete custom_pseq;
            status = FALSE;
            break;
         }
         custom_pseq->add_event((Time_ms)0.0, &pulse1);
         custom_pseq->add_event((Time_ms)TE1/2, &pulse2);
         custom_pseq->add_event(new Sample((Time_ms)TE1,
                                &Phantom::_save_sample, 
                                (void *)phantom));
         custom_pseq->add_event((Time_ms)(TE2+TE1)/2, &pulse3);
         custom_pseq->add_event(new Repeat(TR));
         pseq = (Pulse_Sequence *)custom_pseq;
         break;
      case SCAN_TYPE_DSE_LATE:
         custom_pseq = new Custom_Sequence;
         if (TE2 < TE1){
            cerr << "First echo time must be less than second echo time."
                 << endl;
            delete custom_pseq;
            status = FALSE;
            break;
         }
         custom_pseq->add_event((Time_ms)0.0, &pulse1);
         custom_pseq->add_event((Time_ms)TE1/2, &pulse2);
         custom_pseq->add_event((Time_ms)(TE2+TE1)/2, &pulse3);
         custom_pseq->add_event(new Sample((Time_ms)TE2,
                                &Phantom::_save_sample, 
                                (void *)phantom));
         custom_pseq->add_event(new Repeat(TR));
         pseq = (Pulse_Sequence *)custom_pseq;
         break;
   }
   
   // --- Check for Pulse Sequence creation errors --- //

   if (status == FALSE){
      delete pseq;
      return FALSE;
   }
   paramfile.close();

   // --- Update Pulse Sequence Parameters --- //

   // Update Partial Fourier methods

   Partial_Fourier_Type partial_fourier_method;
   if (scan_percentage != 100.0) {
      partial_fourier_method = PARTIAL_MATRIX;
   } else {
      partial_fourier_method = NO_PARTIAL_FOURIER;
   } 

   if (rectangular_fov != 100.0) {
      partial_fourier_method = RECTANGULAR_FOV;
      scan_percentage = rectangular_fov;
      if (half_scan == YES || partial_echo == YES) {
         cerr << "Can only use one of half scan, partial echo, or "
              << "rectangular FOV at a time." << endl;
         delete pseq; 
         return FALSE;
      }
   } else if (half_scan == YES) {
      partial_fourier_method = HALF_FOURIER;
      if (partial_echo == YES) {
         cerr << "Can only use one of half scan, partial echo, or "
              << "rectangular FOV at a time." << endl;
         delete pseq;
         return FALSE;
      }
   } else if (partial_echo == YES) {
      partial_fourier_method = PARTIAL_ECHO;
   } else {
      partial_fourier_method = NO_PARTIAL_FOURIER;
   }

   // --- Acquisition Geometry --- //

   double voxel_step[3], fov[3];
   int    acquisition_length[3], reconstruction_length[3];

   // Override number of slices defaults
   
   if (args.voxel_dims[SLICE] != 0.0) {
      slice_separation = args.voxel_dims[SLICE];
      slice_thickness  = args.voxel_dims[SLICE];
   }

   if (number_of_slices == 0) {
      number_of_slices = (int)rint((phantom->get_nslices()*
                                    phantom->get_voxel_step(SLICE)) /
                                   slice_separation); 
   }
   if (args.nslices != 0) {
      number_of_slices = args.nslices;
   }

   if (args.slice_thickness != 0.0){
      slice_thickness = args.slice_thickness;
   }

   // Override acquisition length defaults

   if (args.matrix_size != 0) {
      if ((args.matrix_size != 64) && 
          (args.matrix_size != 128) &&
          (args.matrix_size != 256) &&
          (args.matrix_size != 512) &&
          (args.matrix_size != 1024) &&
          (args.matrix_size != 2048) &&
          (args.matrix_size != 4096) &&
          (args.matrix_size != 8192)) {
         cerr << "Matrix size must be one of 64, 128, 256, 512, 1024, 2048, 4096, or 8192." << endl;
         cerr << "Using default: " << scan_matrix << endl;
      } else {
         scan_matrix = args.matrix_size;
      }
   }
 
   acquisition_length[SLICE]  = number_of_slices;
   acquisition_length[ROW]    = scan_matrix;
   acquisition_length[COLUMN] = scan_matrix;

   if (foldover_suppression == TRUE) {
      acquisition_length[foldover_direction] *= 2;
   }

   // Override reconstruction matrix defaults

   if (args.recon_size != 0) {
      if ((args.recon_size != 64) && 
          (args.recon_size != 128) &&
          (args.recon_size != 256) &&
          (args.recon_size != 512) &&
          (args.recon_size != 1024) &&
          (args.recon_size != 2048) &&
          (args.recon_size != 4096) &&
          (args.recon_size != 8192)) {
         cerr << "Reconstruction matrix must be one of 64, 128, 256, 512, 1024, 2048, 4096, or 8192.." 
              << endl;
         cerr << "Using default: " << reconstruction_matrix << endl;
      } else {
         reconstruction_matrix = args.recon_size;
      }
   }

   reconstruction_length[SLICE]  = number_of_slices;
   reconstruction_length[ROW]    = reconstruction_matrix;
   reconstruction_length[COLUMN] = reconstruction_matrix;

   voxel_step[SLICE]  = slice_separation;
   fov[SLICE]  = (number_of_slices-1)*slice_separation + slice_thickness;  

   if ((args.voxel_dims[ROW] != 0.0) && (args.voxel_dims[COLUMN] != 0.0)){
      voxel_step[ROW]    = args.voxel_dims[ROW];
      voxel_step[COLUMN] = args.voxel_dims[COLUMN];
      fov[ROW]    = acquisition_length[ROW]*args.voxel_dims[ROW];
      fov[COLUMN] = acquisition_length[COLUMN]*args.voxel_dims[COLUMN];
   } else {
      voxel_step[ROW]    = field_of_view/acquisition_length[ROW];
      voxel_step[COLUMN] = field_of_view/acquisition_length[COLUMN];
      fov[ROW]    = field_of_view;
      fov[COLUMN] = field_of_view;
   }

   // If old quarter-voxel shift is specified, override the current
   // voxel offset values.

   if (args.oldoffsetFlag) {
      voxel_offset[SLICE]  = (voxel_step[SLICE] -
                              phantom->get_voxel_step(SLICE))/2;
      voxel_offset[ROW]    = (voxel_step[ROW] -
                              phantom->get_voxel_step(ROW))/2;
      voxel_offset[COLUMN] = (voxel_step[COLUMN] -
                              phantom->get_voxel_step(COLUMN))/2;
   }

   pseq->set_slice_thickness(slice_thickness);
   //pseq->set_voxel_step(voxel_step);
   pseq->set_fov(fov);
   pseq->set_voxel_offset(voxel_offset);
   pseq->set_matrix_size(acquisition_length);
   pseq->set_reconstruction_size(reconstruction_length);

   pseq->set_foldover_suppression(foldover_suppression);
   pseq->set_partial_fourier_method(partial_fourier_method, 
                                    scan_percentage/100.0);
   pseq->set_num_of_averages(number_of_signal_averages);
   pseq->set_image_orientation(orientation);
   pseq->set_image_type(image_type);
   pseq->set_scan_mode(scan_mode);
   pseq->set_water_fat_shift(water_fat_shift);

   switch(scan_technique) {
      case SCAN_TYPE_SE:
      case SCAN_TYPE_IR:
      case SCAN_TYPE_SFLASH:
      case SCAN_TYPE_CEFAST:
      case SCAN_TYPE_FISP:
      case SCAN_TYPE_FLASH:
         scanner.apply(quick_pseq);
         break;
      case SCAN_TYPE_DSE_EARLY:
      case SCAN_TYPE_DSE_LATE:
         scanner.apply(custom_pseq);
         break;
   }

   return TRUE;   
}

