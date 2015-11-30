#include <string.h>
#include "mincicv.h"

//==========================================================================
// MINCICV.CXX
// MINC icv class MINC_ICV.
// Inherits from:  
// Base class to:  MINC_File
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

//==========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/mincicv.cxx,v 1.2 2004-08-10 15:23:50 bert Exp $
// $Log: mincicv.cxx,v $
// Revision 1.2  2004-08-10 15:23:50  bert
// Comment out obsolete NC_xxx types
//
// Revision 1.1  2003/05/30 16:43:09  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.5  1996/05/29  16:27:52  rkwan
// Release 2.5
//
// Revision 1.3  1995/12/22  20:17:28  rkwan
// Updated for float type.
//
// Revision 1.2  1995/12/11  14:16:27  rkwan
// Updated for fast_iso_model.
//
//==========================================================================

//--------------------------------------------------------------------------
// MINC_ICV constructors
//--------------------------------------------------------------------------

MINC_ICV::MINC_ICV() {

   int oldncopts =get_ncopts();
   set_ncopts(0);
   if ((_icvid = miicv_create()) == MI_ERROR){
      cerr << "MINC_ICV::MINC_ICV error creating icv." << endl;
   }
   set_ncopts(oldncopts);

}

//--------------------------------------------------------------------------
// MINC_ICV destructors
//--------------------------------------------------------------------------

MINC_ICV::~MINC_ICV() {

#ifdef DEBUG
   assert(_icvid != MI_ERROR);
#endif

   (void)miicv_free(_icvid);
   _icvid = MI_ERROR;

}

//--------------------------------------------------------------------------
// MINC_ICV::get_icv_type
// Sets type_name to a string containing the name of the icv data type.
//--------------------------------------------------------------------------

const char *MINC_ICV::get_icv_type_name(void) {
   static char *type[] = {"byte", "char", "short", "long", "float",
                          "double", "private", "unspecified"};
   char *return_value;
   int  nctype;

   // Get the data type
   if (miicv_inqint(_icvid, MI_ICV_TYPE, &nctype) == MI_ERROR) {
      return (const char *)NULL;
   }
              
   // Set the data type string accordingly
   switch(nctype){
      case NC_BYTE:
         return_value = type[0];
         break;
      case NC_CHAR:
         return_value = type[1];
         break;
      case NC_SHORT:
         return_value = type[2];
         break;
      case NC_LONG:
         return_value = type[3];
         break;
      case NC_FLOAT:
         return_value = type[4];
         break;
      case NC_DOUBLE:
         return_value = type[5];
         break;
#if 0                           // Obsolete
      case NC_BITFIELD:
      case NC_STRING:
      case NC_IARRAY:
      case NC_DIMENSION:
      case NC_VARIABLE:
      case NC_ATTRIBUTE:
      case NC_UNSPECIFIED:
         return_value = type[6];
         break;
#endif /* 0 */
      default:
         return_value = type[7];
         break;
   }

   return return_value;
}

//--------------------------------------------------------------------------
// MINC_ICV::set_default_byte_icv
// Sets up the default icv for a volume with byte data type.
//--------------------------------------------------------------------------

void MINC_ICV::set_default_byte_icv(void){

#ifdef DEBUG
   assert(_icvid != MI_ERROR);
#endif

   // Set desired type
   (void)miicv_setint(_icvid, MI_ICV_TYPE, NC_BYTE);
   (void)miicv_setstr(_icvid, MI_ICV_SIGN, MI_UNSIGNED);

   // Set range of values
   (void)miicv_setint(_icvid, MI_ICV_VALID_MIN, 0);
   (void)miicv_setint(_icvid, MI_ICV_VALID_MAX, 255);

   // Do normalization so that all pixels are on same scale
   (void)miicv_setint(_icvid, MI_ICV_DO_NORM, TRUE);

   // Make sure that any out of range values are mapped to lowest value
   // of type (for input only)
   (void) miicv_setint(_icvid, MI_ICV_DO_FILLVALUE, TRUE);

   // Ensure that images have X, Y, and Z dimensions in the positive
   // direction, giving patient left on left and for drawing from
   // bottom up.  
   (void) miicv_setint(_icvid, MI_ICV_DO_DIM_CONV, TRUE);
   (void) miicv_setint(_icvid, MI_ICV_XDIM_DIR, MI_ICV_POSITIVE);
   (void) miicv_setint(_icvid, MI_ICV_YDIM_DIR, MI_ICV_POSITIVE);
   (void) miicv_setint(_icvid, MI_ICV_ZDIM_DIR, MI_ICV_POSITIVE);

}

//--------------------------------------------------------------------------
// MINC_ICV::set_default_short_icv
// Sets up the default icv for a volume with short data type. 
//--------------------------------------------------------------------------

void MINC_ICV::set_default_short_icv(void){

#ifdef DEBUG
   assert(_icvid != MI_ERROR);
#endif
 
   // Set desired type
   (void)miicv_setint(_icvid, MI_ICV_TYPE, NC_SHORT);
   (void)miicv_setstr(_icvid, MI_ICV_SIGN, MI_SIGNED);

   // Set range of values
   (void)miicv_setint(_icvid, MI_ICV_VALID_MIN, 0);
   (void)miicv_setint(_icvid, MI_ICV_VALID_MAX, 4095);

   // Do normalization so that all pixels are on same scale
   (void)miicv_setint(_icvid, MI_ICV_DO_NORM, TRUE);

   // Make sure that any out of range values are mapped to lowest value
   // of type (for input only)
   (void) miicv_setint(_icvid, MI_ICV_DO_FILLVALUE, TRUE);

   // Ensure that images have X, Y, and Z dimensions in the positive
   // direction, giving patient left on left and for drawing from
   // bottom up.
   (void) miicv_setint(_icvid, MI_ICV_DO_DIM_CONV, TRUE);
   (void) miicv_setint(_icvid, MI_ICV_XDIM_DIR, MI_ICV_POSITIVE);
   (void) miicv_setint(_icvid, MI_ICV_YDIM_DIR, MI_ICV_POSITIVE);
   (void) miicv_setint(_icvid, MI_ICV_ZDIM_DIR, MI_ICV_POSITIVE);

}

//--------------------------------------------------------------------------
// MINC_ICV::set_default_float_icv
// Sets up the default icv for a volume with float data type. 
//--------------------------------------------------------------------------

void MINC_ICV::set_default_float_icv(void){

#ifdef DEBUG
   assert(_icvid != MI_ERROR);
#endif

   // Set desired type
   (void)miicv_setint(_icvid, MI_ICV_TYPE, NC_FLOAT);
   (void)miicv_setstr(_icvid, MI_ICV_SIGN, MI_SIGNED);

   // Do normalization so that all pixels are on same scale
   (void)miicv_setint(_icvid, MI_ICV_DO_NORM, TRUE);

   // Make sure that any out of range values are mapped to lowest value
   // of type (for input only)
   (void) miicv_setint(_icvid, MI_ICV_DO_FILLVALUE, TRUE);

   // Ensure that images have X, Y, and Z dimensions in the positive
   // direction, giving patient left on left and for drawing from
   // bottom up.
   (void) miicv_setint(_icvid, MI_ICV_DO_DIM_CONV, TRUE);
   (void) miicv_setint(_icvid, MI_ICV_XDIM_DIR, MI_ICV_POSITIVE);
   (void) miicv_setint(_icvid, MI_ICV_YDIM_DIR, MI_ICV_POSITIVE);
   (void) miicv_setint(_icvid, MI_ICV_ZDIM_DIR, MI_ICV_POSITIVE);

}
