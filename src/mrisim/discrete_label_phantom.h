#ifndef __DISCRETE_LABEL_PHANTOM_H
#define __DISCRETE_LABEL_PHANTOM_H

//==========================================================================
// DISCRETE_LABEL_PHANTOM.H
// Discrete_Label_Phantom class.
// Inherits from:  Phantom
// Base class to:  Discrete_Phantom, Discrete_RF_Phantom
//
// R. Kwan
// (C) Copyright 1996 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/discrete_label_phantom.h,v 1.1 2003-05-30 16:43:10 bert Exp $
 * $Log: discrete_label_phantom.h,v $
 * Revision 1.1  2003-05-30 16:43:10  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:50:00  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:08:17  rkwan
 * Release 2.5
 *
 *========================================================================*/

#include "phantom.h"
#include <minc/imincfile.h>
#include <minc/omincfile.h>
#include <minc/mrilabel.h>

//---------------------------------------------------------------------------
// Discrete_Label_Phantom class
// Describes a discrete labelled volume MRI phantom.
//---------------------------------------------------------------------------

class Discrete_Label_Phantom : virtual public Phantom {
   public:
      Discrete_Label_Phantom(unsigned int n_tissue_classes); 

      virtual ~Discrete_Label_Phantom();

      // --- Tissue information --- //
      int  open_discrete_label_file(const char *path);
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
      void _load_label_slice(int slice_num, MRI_Label& label_slice);

      // --- Internal data structures --- //
      I_MINC_File _tissue_label_file;

};

//---------------------------------------------------------------------------
// Inline member functions
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::is_same_slice_size_as
// Returns TRUE if the matrix has the same row and column dimensions as
// the phantom.
//---------------------------------------------------------------------------

inline
int Discrete_Label_Phantom::is_same_slice_size_as(const MRI_Matrix& mat) const {
   int same_size;

   if ((mat.get_nrows() == _tissue_label_file.get_nrows()) &&
       (mat.get_ncols() == _tissue_label_file.get_ncols())){
      same_size = TRUE;
   } else {
      same_size = FALSE;
   }

   return same_size;
}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::get_nrows
// Return the length of the ROW dimension of the phantom.
//---------------------------------------------------------------------------

inline
int Discrete_Label_Phantom::get_nrows(void) const {
   return _tissue_label_file.get_nrows();
}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::get_ncols
// Return the length of the COLUMN dimension of the phantom.
//---------------------------------------------------------------------------

inline
int Discrete_Label_Phantom::get_ncols(void) const {
   return _tissue_label_file.get_ncols();
}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::get_nslices
// Return the length of the SLICE dimension of the phantom.
//---------------------------------------------------------------------------

inline
int Discrete_Label_Phantom::get_nslices(void) const {
   return _tissue_label_file.get_nslices();
}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::get_volume_dimensions
// Returns the dimension lengths of the phantom.
//---------------------------------------------------------------------------

inline
void Discrete_Label_Phantom::get_volume_dimensions(int length[]) const {

   length[SLICE]  = _tissue_label_file.get_nslices();
   length[ROW]    = _tissue_label_file.get_nrows();
   length[COLUMN] = _tissue_label_file.get_ncols();

}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::get_volume_info
// Gets information about the labelled phantom volume.
//---------------------------------------------------------------------------

inline
void Discrete_Label_Phantom::get_volume_info(Volume_Info &volume_info) const {
   _tissue_label_file.get_volume_info(volume_info);
}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::get_voxel_step
// Returns the step size of the phantom dimension n.
//---------------------------------------------------------------------------

inline
double Discrete_Label_Phantom::get_voxel_step(int n) const {
   return _tissue_label_file.get_step(n);
}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::get_voxel_start
// Returns the start value of the phantom dimension n.
//---------------------------------------------------------------------------

inline
double Discrete_Label_Phantom::get_voxel_start(int n) const {
   return _tissue_label_file.get_start(n);
}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::display_volume_info
// Displays information about the MINC label volume.
//---------------------------------------------------------------------------

inline
void Discrete_Label_Phantom::display_volume_info(ostream& stream) const {
   stream << "Phantom Labelled Volume Information:" << endl;
   stream << "------------------------------------" << endl << endl;

   _tissue_label_file.display_volume_info(stream);
}

//---------------------------------------------------------------------------
// Discrete_Label_Phantom::set_output_volume_info
// Copies and updates volume information for an output file.
//---------------------------------------------------------------------------

inline
void Discrete_Label_Phantom::set_output_volume_info(O_MINC_File& output,
                                              const Volume_Info &vol_info,
                                              const char *argstring) const {
   output.set_volume_info(_tissue_label_file.get_cdfid(), vol_info,
                          argstring);
}

#endif

