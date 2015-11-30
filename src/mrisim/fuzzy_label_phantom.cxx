//===========================================================================
// FUZZY_LABEL_PHANTOM.CXX
//
// R.Kwan
// April 25, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/fuzzy_label_phantom.cxx,v 1.1 2003-05-30 16:43:10 bert Exp $
// $Log: fuzzy_label_phantom.cxx,v $
// Revision 1.1  2003-05-30 16:43:10  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  15:50:52  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:10:08  rkwan
// Release 2.5
//
//===========================================================================

#include "fuzzy_label_phantom.h"

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom constructor
//---------------------------------------------------------------------------

Fuzzy_Label_Phantom::Fuzzy_Label_Phantom(unsigned int n_tissue_classes) :
   Phantom(n_tissue_classes) {

   _tissue_label_file = new I_MINC_File[n_tissue_classes];

}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom destructor
//---------------------------------------------------------------------------

Fuzzy_Label_Phantom::~Fuzzy_Label_Phantom() {

   this->close_label_files();
   delete[] _tissue_label_file;

}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::open_fuzzy_label_file
//---------------------------------------------------------------------------

int Fuzzy_Label_Phantom::open_fuzzy_label_file(Tissue_Label tissue_label,
                                         const char *path) {

   unsigned int tissue_index = get_tissue_index(tissue_label);
   I_MINC_File &label_file   = _tissue_label_file[tissue_index];

   // Open the MINC file
   label_file.open(path);

   // Set up the ICV so that scaling is done to the maximum tissue id
   // label in the phantom
   label_file.set_icv_property(MI_ICV_TYPE, (int)NC_FLOAT);
   label_file.set_icv_property(MI_ICV_DO_NORM, TRUE);
   label_file.set_icv_property(MI_ICV_DO_FILLVALUE, TRUE);
   label_file.set_icv_property(MI_ICV_DO_DIM_CONV, TRUE);
   label_file.set_icv_property(MI_ICV_XDIM_DIR, MI_ICV_POSITIVE);
   label_file.set_icv_property(MI_ICV_YDIM_DIR, MI_ICV_POSITIVE);
   label_file.set_icv_property(MI_ICV_ZDIM_DIR, MI_ICV_POSITIVE);

   // Attach the ICV 
   label_file.attach_icv();

   // Check label file size for consistency

   int consistent = ((!_tissue_label_file[0].is_good()) ? TRUE : 
                     ((label_file.get_nrows() == get_nrows()) &&
                      (label_file.get_ncols() == get_ncols()) &&
                      (label_file.get_nslices() == get_nslices())));

   if (!consistent) {
      cerr << endl << "Fuzzy volume size is inconsistent with phantom: "
           << _tissue_label_file[0].get_filename() << "."
           << endl; 
   }
      
   return (label_file.is_good() && consistent);

}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::close_label_files
//---------------------------------------------------------------------------

void Fuzzy_Label_Phantom::close_label_files(void){

   unsigned int itissue;
   I_MINC_File  *label_file;

   for (itissue=0; itissue<get_num_tissues(); itissue++){
      label_file = &(_tissue_label_file[itissue]);

      if (label_file->is_open()){
         label_file->close();
      }
   }

}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::_load_label_slice
// Load a slice of the labelled volume into memory from a MINC file.
//---------------------------------------------------------------------------

void Fuzzy_Label_Phantom::_load_label_slice(int slice_num,
                                      Tissue_Label tissue_label,
                                      Real_Slice& label_slice) {

#ifdef DEBUG
   assert(this->is_same_slice_size_as(label_slice));
#endif

   unsigned int tissue_index = Phantom::get_tissue_index(tissue_label); 
   int r = _tissue_label_file[tissue_index].load_slice(slice_num, (void *)label_slice);
   if (r < 0) {
     cerr << endl << "Failed to load slice #" << slice_num;
     cerr << " for tissue label " << (int)tissue_label << "." << endl;
     cerr << "You may have omitted the input file for that tissue type." << endl;
     exit(EXIT_FAILURE);
   }

}
