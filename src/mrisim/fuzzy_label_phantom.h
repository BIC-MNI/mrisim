#ifndef __FUZZY_LABEL_PHANTOM_H
#define __FUZZY_LABEL_PHANTOM_H

//==========================================================================
// FUZZY_LABEL_PHANTOM.H
// Fuzzy_Label_Phantom class.
// Inherits from:  Phantom
// Base class to:  Fuzzy_Phantom, Fuzzy_RF_Phantom
//
// R. Kwan
// (C) Copyright 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/fuzzy_label_phantom.h,v 1.1 2003-05-30 16:43:10 bert Exp $
 * $Log: fuzzy_label_phantom.h,v $
 * Revision 1.1  2003-05-30 16:43:10  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:50:46  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:10:00  rkwan
 * Release 2.5
 *
 *========================================================================*/

#include "phantom.h"
#include <minc/imincfile.h>
#include <minc/mrimatrix.h>

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom class
// Describes a fuzzy labelled MRI phantom.
//---------------------------------------------------------------------------

class Fuzzy_Label_Phantom : virtual public Phantom {
   public:
      Fuzzy_Label_Phantom(unsigned int n_tissue_classes);

      virtual ~Fuzzy_Label_Phantom();

      // --- Tissue Information --- //
      int  open_fuzzy_label_file(Tissue_Label tissue_label,
                                 const char *path);
      void close_label_files(void);

      // --- Volume convenience functions --- //
      inline int    is_same_slice_size_as(const MRI_Matrix& mat) const;
      inline int    get_nrows(void) const;
      inline int    get_ncols(void) const;
      inline int    get_nslices(void) const;
      inline void   get_volume_dimensions(int length[]) const;
      inline void   get_volume_info(Volume_Info &volume_info) const;
      inline double get_voxel_step(int n) const;
      inline double get_voxel_start(int n) const;
      inline void   display_volume_info(ostream& stream) const;
      inline void   set_output_volume_info(O_MINC_File& output,
                                        const Volume_Info& vol_info,
                                        const char *argstring = NULL) const;
 
   protected:

      // --- Internal member functions --- //
      void _load_label_slice(int slice_num, 
                             Tissue_Label tissue_label,
                             Real_Slice& label_slice);

      // --- Internal data structures --- //
      I_MINC_File *_tissue_label_file;

};

//---------------------------------------------------------------------------
// Inline member functions
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::is_same_slice_size_as
// Returns TRUE if the matrix has the same row and column dimensions as
// the phantom.
//---------------------------------------------------------------------------

inline
int Fuzzy_Label_Phantom::is_same_slice_size_as(const MRI_Matrix& mat) const {
   int same_size;

   if ((mat.get_nrows() == _tissue_label_file[0].get_nrows()) &&
       (mat.get_ncols() == _tissue_label_file[0].get_ncols())){
      same_size = TRUE;
   } else {
      same_size = FALSE;
   }

   return same_size;
}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::get_nrows
// Return the length of the ROW dimension of the phantom.
//---------------------------------------------------------------------------

inline
int Fuzzy_Label_Phantom::get_nrows(void) const {
   return _tissue_label_file[0].get_nrows();
}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::get_ncols
// Return the length of the COLUMN dimension of the phantom.
//---------------------------------------------------------------------------

inline
int Fuzzy_Label_Phantom::get_ncols(void) const {
   return _tissue_label_file[0].get_ncols();
}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::get_nslices
// Return the length of the SLICE dimension of the phantom.
//---------------------------------------------------------------------------

inline
int Fuzzy_Label_Phantom::get_nslices(void) const {
   return _tissue_label_file[0].get_nslices();
}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::get_volume_dimensions
// Returns the dimension lengths of the phantom.
//---------------------------------------------------------------------------

inline
void Fuzzy_Label_Phantom::get_volume_dimensions(int length[]) const {

   length[SLICE]  = _tissue_label_file[0].get_nslices();
   length[ROW]    = _tissue_label_file[0].get_nrows();
   length[COLUMN] = _tissue_label_file[0].get_ncols();

}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::get_volume_info
// Gets information about the labelled phantom volume.
//---------------------------------------------------------------------------

inline
void Fuzzy_Label_Phantom::get_volume_info(Volume_Info &volume_info) const {
   _tissue_label_file[0].get_volume_info(volume_info);
}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::get_voxel_step
// Returns the step size of the phantom dimension n.
//---------------------------------------------------------------------------

inline
double Fuzzy_Label_Phantom::get_voxel_step(int n) const {
   return _tissue_label_file[0].get_step(n);
}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::get_voxel_start
// Returns the start value of the phantom dimension n.
//---------------------------------------------------------------------------

inline
double Fuzzy_Label_Phantom::get_voxel_start(int n) const {
   return _tissue_label_file[0].get_start(n);
}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::display_volume_info 
// Displays information about the MINC label volume.
//---------------------------------------------------------------------------

inline
void Fuzzy_Label_Phantom::display_volume_info(ostream& stream) const {
   stream << "Phantom Labelled Volume Information:" << endl;
   stream << "------------------------------------" << endl << endl;
   _tissue_label_file[0].display_volume_info(stream);
}

//---------------------------------------------------------------------------
// Fuzzy_Label_Phantom::set_output_volume_info
// Copies and updates volume information for an output file.
//---------------------------------------------------------------------------

inline
void Fuzzy_Label_Phantom::set_output_volume_info(O_MINC_File& output,
                                           const Volume_Info& vol_info,
                                           const char *argstring) const {

   output.set_volume_info(_tissue_label_file[0].get_cdfid(), vol_info,
                          argstring);

}

#endif
