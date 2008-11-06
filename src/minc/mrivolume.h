#ifndef __MRIVOLUME_H
#define __MRIVOLUME_H

//===========================================================================
// MRIVOLUME.H
//
// R.Kwan
// August 31, 1995
//
// (C) Copyright 1995 by R.Kwan
//===========================================================================

/*===========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/mrivolume.h,v 1.3 2008-11-06 10:58:23 rotor Exp $
 * $Log: mrivolume.h,v $
 * Revision 1.3  2008-11-06 10:58:23  rotor
 *  * fixed includes for iostream and friends
 *  * updated for new release (1.0.2)
 *
 * Revision 1.2  2004/08/10 15:35:37  bert
 * Add 'class' keyword to friend declarations
 *
 * Revision 1.1  2003/05/30 16:43:09  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:28:55  rkwan
 * Release 2.5
 *
 * Revision 1.5  1996/01/08  15:33:39  rkwan
 * Added MRI_Float_Volume.
 *
 * Revision 1.4  1995/12/11  15:14:33  rkwan
 * Fix RCS header bug.
 *
 *=========================================================================*/

#include "mrimatrix.h"

#include <iostream>
#include "fourn.h"

#ifdef DEBUG
#include <assert.h>
#endif

//===========================================================================
// Base MRI_Volume class
// Abstract class which defines common interface to matrices.
//===========================================================================

class MRI_Volume {
   public:
      MRI_Volume();
      MRI_Volume(unsigned int nrows, unsigned int ncols, unsigned int nslices);
 
      virtual ~MRI_Volume();

      unsigned int get_nrows(void) const {return _nrows;}
      unsigned int get_ncols(void) const {return _ncols;}
      unsigned int get_nslices(void) const {return _nslices;}
      unsigned int get_nelements(void) const {return _nrows*_ncols*_nslices;}
      unsigned int get_slicesize(void) const {return _nrows*_ncols;}

      virtual size_t size_in_bytes(void) const = 0;
      virtual size_t element_size_in_bytes(void) const = 0;

      void set_nrows(unsigned int nrows) {_nrows = nrows;}
      void set_ncols(unsigned int ncols) {_ncols = ncols;}
      void set_nslices(unsigned int nslices) {_nslices = nslices;}

      void reshape(unsigned int nrows, unsigned int ncols, 
                   unsigned int nslices){
         if ((nrows*ncols*nslices) == (_nrows*_ncols*_nslices)){
            _nrows = nrows; _ncols = ncols; _nslices = nslices; }
      }

      int is_same_size_as(MRI_Volume& vol) const;
      int slice_size_is_same_as(MRI_Matrix& mat) const;
      int is_power_of_two(void) const;

      unsigned int get_offset(unsigned int row, unsigned int col, 
                              unsigned int slice){
         return (slice*_nrows + row)*_ncols + col; }

      virtual operator void *() = 0;

   protected:
  
      unsigned int _nrows;
      unsigned int _ncols;
      unsigned int _nslices;

      virtual void _allocate(unsigned int nelements) = 0;
      virtual void _deallocate(void) = 0;

};

//---------------------------------------------------------------------------
// Forward references of type-dependent matrix classes.
//---------------------------------------------------------------------------

class MRI_Byte_Volume;
class MRI_Short_Volume;
class MRI_Float_Volume;
class MRI_Double_Volume;
class MRI_Complex_Volume;

//===========================================================================
// MRI_Byte_Volume
//===========================================================================

class MRI_Byte_Volume : public MRI_Volume {
   public:
      MRI_Byte_Volume();
      MRI_Byte_Volume(unsigned int nrows, unsigned int ncols, 
                    unsigned int nslices, unsigned char fill_value=0);
      MRI_Byte_Volume(const MRI_Byte_Volume& vol);
      MRI_Byte_Volume(const MRI_Short_Volume& vol);
      MRI_Byte_Volume(const MRI_Float_Volume& vol);
      MRI_Byte_Volume(const MRI_Double_Volume& vol);

      virtual ~MRI_Byte_Volume();
 
      virtual size_t size_in_bytes(void) const {
         return get_nelements()*sizeof(unsigned char); }
      virtual size_t element_size_in_bytes(void) const {
         return sizeof(unsigned char); }

      void fill_with(unsigned char fill_value);
      void zeros(void) { fill_with(0); }

      MRI_Byte_Volume& operator=(const MRI_Byte_Volume& vol);
      MRI_Byte_Volume& operator*=(double a); 

      int operator != (const MRI_Byte_Volume& vol) const;
      int operator == (const MRI_Byte_Volume& vol) const {
         return !(operator != (vol)); }

      unsigned char maximum(void) const;
      unsigned char minimum(void) const;
      double sum(void) const;
      double norm(void) const;
      double mean(void) const;
      double std(void) const;

      virtual operator void *() {return (void *)_volume;}
      operator unsigned char *() {return _volume;}
      unsigned char *operator[] (int slice){
         return &(_volume[slice*_nrows*_ncols]); }

      unsigned char& operator() (unsigned int row, unsigned int col,
                                 unsigned int slice){
         return _volume[(slice*_nrows + row)*_ncols + col]; }
      unsigned char& operator() (unsigned int row, unsigned int col,
                                 unsigned int slice) const {
         return _volume[(slice*_nrows + row)*_ncols + col]; }
      
      MRI_Byte_Matrix operator() (unsigned int slice) const;
      MRI_Byte_Volume operator() (unsigned int row1, unsigned int row2,
                                unsigned int col1, unsigned int col2, 
                                unsigned int slice1, unsigned int slice2) const;

      void display(ostream& stream);

   protected:
      friend class MRI_Short_Volume;
      friend class MRI_Float_Volume;
      friend class MRI_Double_Volume;
  
      unsigned char *_volume;
  
      virtual void _allocate(unsigned int nelements) {
         _volume = new unsigned char[nelements]; }
      virtual void _deallocate(void) {
         delete[] _volume; }
      
};

//===========================================================================
// MRI_Short_Volume
//===========================================================================

class MRI_Short_Volume : public MRI_Volume {
   public:
      MRI_Short_Volume();
      MRI_Short_Volume(unsigned int nrows, unsigned int ncols,
                     unsigned int nslices, short fill_value=0);
      MRI_Short_Volume(const MRI_Byte_Volume& vol);
      MRI_Short_Volume(const MRI_Short_Volume& vol);
      MRI_Short_Volume(const MRI_Float_Volume& vol);
      MRI_Short_Volume(const MRI_Double_Volume& vol);

      virtual ~MRI_Short_Volume();

      virtual size_t size_in_bytes(void) const {
         return get_nelements()*sizeof(short); }
      virtual size_t element_size_in_bytes(void) const {
         return sizeof(short); }

      void fill_with(short fill_value);
      void zeros(void) { fill_with(0); }

      MRI_Short_Volume& operator=(const MRI_Short_Volume& vol);
      MRI_Short_Volume& operator*=(double a);

      int operator != (const MRI_Short_Volume& vol) const;
      int operator == (const MRI_Short_Volume& vol) const {
         return !(operator != (vol)); }

      short maximum(void) const;
      short minimum(void) const;
      double sum(void) const;
      double norm(void) const;
      double mean(void) const;
      double std(void) const;

      virtual operator void *() {return (void *)_volume;}
      operator short *() {return _volume;}
      short *operator[] (int slice){
         return &(_volume[slice*_nrows*_ncols]); }

      short& operator() (unsigned int row, unsigned int col,
                                 unsigned int slice){
         return _volume[(slice*_nrows + row)*_ncols + col]; }
      short& operator() (unsigned int row, unsigned int col,
                                 unsigned int slice) const {
         return _volume[(slice*_nrows + row)*_ncols + col]; }

      MRI_Short_Matrix operator() (unsigned int slice) const;
      MRI_Short_Volume operator() (unsigned int row1, unsigned int row2,
                                unsigned int col1, unsigned int col2,
                                unsigned int slice1, unsigned int slice2) const;

      void display(ostream& stream);

   protected:
      friend class MRI_Byte_Volume;
      friend class MRI_Float_Volume;
      friend class MRI_Double_Volume;

      short *_volume;

      virtual void _allocate(unsigned int nelements) {
         _volume = new short[nelements]; }
      virtual void _deallocate(void) {
         delete[] _volume; }

};

//===========================================================================
// MRI_Float_Volume
//===========================================================================

class MRI_Float_Volume : public MRI_Volume {
   public:
      MRI_Float_Volume();
      MRI_Float_Volume(unsigned int nrows, unsigned int ncols, 
                      unsigned int nslices, float fill = 0.0);
      MRI_Float_Volume(const MRI_Float_Volume& vol);
      MRI_Float_Volume(const MRI_Byte_Volume& vol);
      MRI_Float_Volume(const MRI_Short_Volume& vol);
      MRI_Float_Volume(const MRI_Double_Volume& vol);

      virtual ~MRI_Float_Volume();

      virtual size_t size_in_bytes(void) const {
         return get_nelements()*sizeof(float); }
      virtual size_t element_size_in_bytes(void) const {
         return sizeof(float); }

      void fill_with(float fill_value);
      void zeros(void) { fill_with(0); }

      MRI_Float_Volume& operator=(const MRI_Float_Volume& vol);
      MRI_Float_Volume& operator*=(double a);

      int operator != (const MRI_Float_Volume& vol) const;
      int operator == (const MRI_Float_Volume& vol) const {
         return !(operator != (vol)); }

      float maximum(void) const;
      float minimum(void) const;
      double sum(void) const;
      double norm(void) const;
      double mean(void) const;
      double std(void) const;

      MRI_Complex_Volume FFT2(void) const;
      MRI_Complex_Volume iFFT2(void) const;
      MRI_Complex_Volume FFT3(void) const;
      MRI_Complex_Volume iFFT3(void) const;
      void fftshift2(void);
      void fftshift3(void);

      virtual operator void *() {return (void *)_volume;}
      operator float *() {return _volume;}
      float *operator [] (int slice){
         return &(_volume[slice*_nrows*_ncols]); }

      float& operator() (unsigned int row, unsigned int col, 
                         unsigned int slice){
         return _volume[(slice*_nrows + row)*_ncols + col]; }

      float& operator() (unsigned int row, unsigned int col, 
                         unsigned int slice) const {
         return _volume[(slice*_nrows + row)*_ncols + col]; }

      MRI_Float_Matrix operator() (unsigned int slice) const;         
      MRI_Float_Volume operator() (unsigned int row1, unsigned int row2, 
                                   unsigned int row3,
                                   unsigned int col1, unsigned int col2, 
                                   unsigned int col3) const;

      void display(ostream& stream);

   protected:
      friend class MRI_Byte_Volume;
      friend class MRI_Short_Volume;
      friend class MRI_Double_Volume;
      friend class MRI_Complex_Volume;

      float *_volume;
  
      virtual void _allocate(unsigned int nelements) {
         _volume = new float[nelements]; }
      virtual void _deallocate(void) {
         delete[] _volume; }
};

//===========================================================================
// MRI_Double_Volume
//===========================================================================

class MRI_Double_Volume : public MRI_Volume {
   public:
      MRI_Double_Volume();
      MRI_Double_Volume(unsigned int nrows, unsigned int ncols, 
                      unsigned int nslices, double fill = 0.0);
      MRI_Double_Volume(const MRI_Double_Volume& vol);
      MRI_Double_Volume(const MRI_Byte_Volume& vol);
      MRI_Double_Volume(const MRI_Short_Volume& vol);

      virtual ~MRI_Double_Volume();

      virtual size_t size_in_bytes(void) const {
         return get_nelements()*sizeof(double); }
      virtual size_t element_size_in_bytes(void) const {
         return sizeof(double); }

      void fill_with(double fill_value);
      void zeros(void) { fill_with(0); }

      MRI_Double_Volume& operator=(const MRI_Double_Volume& vol);
      MRI_Double_Volume& operator*=(double a);

      int operator != (const MRI_Double_Volume& vol) const;
      int operator == (const MRI_Double_Volume& vol) const {
         return !(operator != (vol)); }

      double maximum(void) const;
      double minimum(void) const;
      double sum(void) const;
      double norm(void) const;
      double mean(void) const;
      double std(void) const;

      MRI_Complex_Volume FFT2(void) const;
      MRI_Complex_Volume iFFT2(void) const;
      MRI_Complex_Volume FFT3(void) const;
      MRI_Complex_Volume iFFT3(void) const;
      void fftshift2(void);
      void fftshift3(void);

      virtual operator void *() {return (void *)_volume;}
      operator double *() {return _volume;}
      double *operator [] (int slice){
         return &(_volume[slice*_nrows*_ncols]); }

      double& operator() (unsigned int row, unsigned int col, 
                          unsigned int slice){
         return _volume[(slice*_nrows + row)*_ncols + col]; }

      double& operator() (unsigned int row, unsigned int col, 
                          unsigned int slice) const {
         return _volume[(slice*_nrows + row)*_ncols + col]; }

      MRI_Double_Matrix operator() (unsigned int slice) const;         
      MRI_Double_Volume operator() (unsigned int row1, unsigned int row2, 
                                    unsigned int row3,
                                    unsigned int col1, unsigned int col2, 
                                    unsigned int col3) const;

      void display(ostream& stream);

   protected:
      friend class MRI_Byte_Volume;
      friend class MRI_Short_Volume;
      friend class MRI_Float_Volume;
      friend class MRI_Complex_Volume;

      double *_volume;
  
      virtual void _allocate(unsigned int nelements) {
         _volume = new double[nelements]; }
      virtual void _deallocate(void) {
         delete[] _volume; }
};

//===========================================================================
// MRI_Complex_Volume
//===========================================================================

class MRI_Complex_Volume : public MRI_Volume {
   public:
      MRI_Complex_Volume();
      MRI_Complex_Volume(unsigned int nrows, unsigned int ncols, 
                         unsigned int nslices,
                         double fill = 0.0);
      MRI_Complex_Volume(const MRI_Complex_Volume& vol);
      MRI_Complex_Volume(const MRI_Float_Volume& vol);
      MRI_Complex_Volume(const MRI_Double_Volume& vol);

      virtual ~MRI_Complex_Volume();

      virtual size_t size_in_bytes(void) const {
         return 2*get_nelements()*sizeof(double); }
      virtual size_t element_size_in_bytes(void) const {
         return 2*sizeof(double); }
     
      void fill_with(double fill_value);
      void fill_real_with(double fill_value);
      void fill_imag_with(double fill_value);
      void zeros(void) { fill_with(0); }

      MRI_Double_Volume real(void) const;
      MRI_Double_Volume imag(void) const;
      MRI_Double_Volume abs(void)  const;
      MRI_Double_Volume angle(void) const;

      MRI_Complex_Volume& operator=(const MRI_Complex_Volume& vol);
      MRI_Complex_Volume& operator*=(double a);

      void FFT2(void);
      void iFFT2(void);
      void FFT3(void);
      void iFFT3(void);
      void fftshift2(void);
      void fftshift3(void);

      virtual operator void *() {return (void *)_volume;}
      operator double *() {return _volume;}

      double& real(unsigned int row, unsigned int col, unsigned int slice){
         return _volume[2*((slice*_nrows + row)*_ncols + col)]; }
      double& imag(unsigned int row, unsigned int col, unsigned int slice){
         return _volume[2*((slice*_nrows + row)*_ncols + col)+1]; }

      void display(ostream& stream);

   protected:
      
      double *_volume;

      virtual void _allocate(unsigned int nelements) {
         _volume = new double[2*nelements]; }
      virtual void _deallocate(void) {
         delete[] _volume; }


};

#endif
