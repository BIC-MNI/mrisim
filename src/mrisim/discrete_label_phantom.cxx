//===========================================================================
// DISCRETE_LABEL_PHANTOM.CXX
//
// R.Kwan
// April 25, 1996
//
// (C) Copyright 1996 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/mrisim/discrete_label_phantom.cxx,v 1.1 2003-05-30 16:43:10 bert Exp $
// $Log: discrete_label_phantom.cxx,v $
// Revision 1.1  2003-05-30 16:43:10  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  15:50:11  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:08:44  rkwan
// Release 2.5
//
//===========================================================================

#include "discrete_label_phantom.h"

//---------------------------------------------------------------------------
// Discrete_Label_Phantom constructor
//---------------------------------------------------------------------------

Discrete_Label_Phantom::Discrete_Label_Phantom(unsigned int n_tissue_classes) :
   Phantom(n_tissue_classes) {

}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom destructor
//---------------------------------------------------------------------------

Discrete_Label_Phantom::~Discrete_Label_Phantom() {

   this->close_label_files();

}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::open_discrete_label_file
//---------------------------------------------------------------------------

int Discrete_Label_Phantom::open_discrete_label_file(const char *path) {

   // Open the MINC file
   _tissue_label_file.open(path);

   // Set up the ICV so that scaling is done to the maximum tissue id
   // label in the phantom
   _tissue_label_file.set_icv_property(MI_ICV_TYPE, (int)NC_BYTE);
   _tissue_label_file.set_icv_property(MI_ICV_SIGN, MI_UNSIGNED);
   _tissue_label_file.set_icv_property(MI_ICV_VALID_MIN, 0);
   _tissue_label_file.set_icv_property(MI_ICV_VALID_MAX, 255);
   _tissue_label_file.set_icv_property(MI_ICV_DO_NORM, FALSE);
   _tissue_label_file.set_icv_property(MI_ICV_DO_FILLVALUE, TRUE);
   _tissue_label_file.set_icv_property(MI_ICV_DO_DIM_CONV, TRUE);
   _tissue_label_file.set_icv_property(MI_ICV_XDIM_DIR, MI_ICV_POSITIVE);
   _tissue_label_file.set_icv_property(MI_ICV_YDIM_DIR, MI_ICV_POSITIVE);
   _tissue_label_file.set_icv_property(MI_ICV_ZDIM_DIR, MI_ICV_POSITIVE);

   // Attach the ICV 
   _tissue_label_file.attach_icv();

   return _tissue_label_file.is_good();

}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::close_label_files
//---------------------------------------------------------------------------

void Discrete_Label_Phantom::close_label_files(void) {

   if (_tissue_label_file.is_open()){
      _tissue_label_file.close();
   }

}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::_load_label_slice
// Load a slice of the labelled volume into memory from a MINC file.
//---------------------------------------------------------------------------

void Discrete_Label_Phantom::_load_label_slice(int slice_num, 
                                         MRI_Label& label_slice) {

#ifdef DEBUG
   // Check that the label slice is the correct size for the phantom
   assert(this->is_same_slice_size_as(label_slice));
#endif

   // Read in the real image-min and image-max for the slice
   long start[3] = {0, 0, 0};
   start[SLICE]  = slice_num;

   double image_min, image_max;

   _tissue_label_file.get_variable(MIimagemin, start, NC_DOUBLE, NULL,
                                   &image_min);
   _tissue_label_file.get_variable(MIimagemax, start, NC_DOUBLE, NULL,
                                   &image_max);

   // Update the ICV to read labels with the correct real image-min
   // and image-max values

   _tissue_label_file.detach_icv();
   _tissue_label_file.set_icv_property(MI_ICV_VALID_MIN,
                                   (unsigned char)rint(image_min));
   _tissue_label_file.set_icv_property(MI_ICV_VALID_MAX,
                                   (unsigned char)rint(image_max));
   _tissue_label_file.attach_icv();

   // Load the phantom slice into memory

   _tissue_label_file.load_slice(slice_num, (void *)label_slice);

}
