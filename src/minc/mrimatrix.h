#ifndef __MRIMATRIX_H
#define __MRIMATRIX_H

//===========================================================================
// MRIMATRIX.H
// Matrix support classes.
// Inherits from:
// Base class to:  MRI_Image, MRI_Volume
// R.Kwan
// August 31, 1995
//
// (C) Copyright 1995 by R.Kwan
//===========================================================================

/*===========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/minc/mrimatrix.h,v 1.3 2008-11-06 10:58:23 rotor Exp $
 * $Log: mrimatrix.h,v $
 * Revision 1.3  2008-11-06 10:58:23  rotor
 *  * fixed includes for iostream and friends
 *  * updated for new release (1.0.2)
 *
 * Revision 1.2  2004/08/10 15:35:14  bert
 * Add 'class' keyword to friend declarations, use atan2() rather than fatan2()
 *
 * Revision 1.1  2003/05/30 16:43:09  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:48:29  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:24:52  rkwan
 * Release 2.5
 *
 * Revision 3.3  1995/12/22  20:18:07  rkwan
 * Updated for float type.
 *
 * Revision 3.2  1995/12/11  15:13:58  rkwan
 * Fix RCS header bug.
 *
 *=========================================================================*/

#include <math.h>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <float.h>
#include <string.h>
#include "fourn.h"

using namespace std;

#ifdef DEBUG
#include <assert.h>
#include <stdio.h>
#endif

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

//===========================================================================
// Base MRI_Matrix class
// Abstract class which defines common interface to matrices.
//===========================================================================

class MRI_Matrix {
   public:
      MRI_Matrix(unsigned int nrows, unsigned int ncols);
 
      virtual ~MRI_Matrix();

      unsigned int get_nrows(void) const {return _nrows;}
      unsigned int get_ncols(void) const {return _ncols;}
      unsigned int get_nelements(void) const {return _nrows*_ncols;}

      virtual size_t size_in_bytes(void) const = 0;
      virtual size_t element_size_in_bytes(void) const = 0;

      void set_nrows(unsigned int nrows) {_nrows = nrows;}
      void set_ncols(unsigned int ncols) {_ncols = ncols;}

      int reshape(unsigned int nrows, unsigned int ncols) {
         int status;
         if ((nrows*ncols) == (_nrows*_ncols)){
            _nrows = nrows; _ncols = ncols; status = TRUE;
         } else {
            status = FALSE;
         }
         return status;
      }

      int is_same_size_as(const MRI_Matrix& mat) const;
      int is_power_of_two(void) const;

      unsigned int get_offset(unsigned int row, unsigned int col) const {
         return row*_ncols + col;
      }
      unsigned int get_row_stride(void) const {
         return (unsigned int)1; 
      }
      unsigned int get_col_stride(void) const {
         return _ncols;
      }
      unsigned int get_row_length(void) const {
         return _ncols;
      }
      unsigned int get_col_length(void) const {
         return _nrows;
      }

      virtual operator void *() = 0;
      virtual operator const void *() const = 0;

   protected:
 
      MRI_Matrix();
 
      unsigned int _nrows;
      unsigned int _ncols;

      virtual void _allocate(unsigned int nelements) = 0;
      virtual void _deallocate(void) = 0;

};

//---------------------------------------------------------------------------
// Forward references of type-dependent matrix classes.
//---------------------------------------------------------------------------

class MRI_Byte_Matrix;
class MRI_Short_Matrix;
class MRI_Float_Matrix;
class MRI_Double_Matrix;
class MRI_FComplex_Matrix;
class MRI_Complex_Matrix;
class MRI_Label;
class MRI_Image;

//===========================================================================
// MRI_Byte_Matrix
//===========================================================================

class MRI_Byte_Matrix : public MRI_Matrix {
   public:
      MRI_Byte_Matrix(unsigned int nrows, unsigned int ncols, 
                      unsigned char fill=0);
      MRI_Byte_Matrix(const MRI_Byte_Matrix& mat);
      MRI_Byte_Matrix(const MRI_Short_Matrix& mat);
      MRI_Byte_Matrix(const MRI_Float_Matrix& mat);
      MRI_Byte_Matrix(const MRI_Double_Matrix& mat);

      virtual ~MRI_Byte_Matrix();

      virtual size_t size_in_bytes(void) const {
         return get_nelements()*sizeof(unsigned char); }
      virtual size_t element_size_in_bytes(void) const {
         return sizeof(unsigned char); }
 
      void fill_with(unsigned char fill_value);
      void zeros(void) { fill_with(0); }
      void ones(void)  { fill_with(1); }

      MRI_Byte_Matrix& operator=(const MRI_Byte_Matrix& mat);
      MRI_Byte_Matrix& operator*=(double a); 
      MRI_Byte_Matrix& operator*=(const MRI_Byte_Matrix& x);
      MRI_Byte_Matrix& operator+=(const MRI_Byte_Matrix& x);
      MRI_Byte_Matrix& saxpy(double a, const MRI_Byte_Matrix& x,
                             const MRI_Byte_Matrix& y);
      MRI_Byte_Matrix& gaxpy(const MRI_Byte_Matrix& a, 
                             const MRI_Byte_Matrix& x,
                             const MRI_Byte_Matrix& y);

      int operator != (const MRI_Byte_Matrix& mat) const;
      int operator == (const MRI_Byte_Matrix& mat) const {
         return !(operator != (mat)); }

      unsigned char maximum(void) const;
      unsigned char minimum(void) const;
      double sum(void) const;
      double norm(void) const;
      double mean(void) const;
      double std(void) const;

      unsigned char maximum(unsigned int row1, unsigned int row2,
                            unsigned int col1, unsigned int col2) const;
      unsigned char minimum(unsigned int row1, unsigned int row2,
                            unsigned int col1, unsigned int col2) const;
      double sum(unsigned int row1, unsigned int row2,
                 unsigned int col1, unsigned int col2) const;
      double norm(unsigned int row1, unsigned int row2,
                  unsigned int col1, unsigned int col2) const;
      double mean(unsigned int row1, unsigned int row2,
                  unsigned int col1, unsigned int col2) const;
      double std(unsigned int row1, unsigned int row2,
                 unsigned int col1, unsigned int col2) const;

      void set_submatrix(MRI_Byte_Matrix& mat, 
                         unsigned int row, unsigned int col);

      unsigned char *row_ptr(int row){
         return &(_matrix[row*_ncols]); }
      unsigned char *col_ptr(int col){
         return &(_matrix[col]); }

      const unsigned char *row_ptr(int row) const {
         return &(_matrix[row*_ncols]); }
      const unsigned char *col_ptr(int col) const {
         return &(_matrix[col]); }

      virtual operator void *() {
         return (void *)_matrix; }
      virtual operator const void *() const {
         return (const void *)_matrix; }

      operator unsigned char *() {
         return _matrix; }
      operator const unsigned char *() const {
         return _matrix; }

      unsigned char *operator[] (int row){
         return row_ptr(row); }
      const unsigned char *operator[] (int row) const {
         return row_ptr(row); }

      unsigned char& operator() (unsigned int row, unsigned int col){
         return _matrix[row*_ncols+col]; }
      unsigned char  operator() (unsigned int row, unsigned int col) const {
         return _matrix[row*_ncols+col]; }
      MRI_Byte_Matrix operator() (unsigned int row1, unsigned int row2,
                                  unsigned int col1, unsigned int col2) const;

      void display(ostream& stream) const;

   protected:
      MRI_Byte_Matrix();

      unsigned char *_matrix;

      const unsigned char *_matrix_endptr(void) const {
         return _matrix+get_nelements(); }

      virtual void _allocate(unsigned int nelements) {
         _matrix = new unsigned char[nelements]; }
      virtual void _deallocate(void) {
         delete[] _matrix; }
      
   private:
      friend class MRI_Short_Matrix;
      friend class MRI_Float_Matrix;
      friend class MRI_Double_Matrix;
      friend class MRI_FComplex_Matrix;
      friend class MRI_Complex_Matrix;
      friend class MRI_Label;

};

//===========================================================================
// MRI_Short_Matrix
//===========================================================================

class MRI_Short_Matrix : public MRI_Matrix {
   public:
      MRI_Short_Matrix(unsigned int nrows, unsigned int ncols, 
                       short fill=0);
      MRI_Short_Matrix(const MRI_Byte_Matrix& mat);
      MRI_Short_Matrix(const MRI_Short_Matrix& mat);
      MRI_Short_Matrix(const MRI_Float_Matrix& mat);
      MRI_Short_Matrix(const MRI_Double_Matrix& mat);

      virtual ~MRI_Short_Matrix();

      virtual size_t size_in_bytes(void) const {
         return get_nelements()*sizeof(short); }
      virtual size_t element_size_in_bytes(void) const {
         return sizeof(short); }
 
      void fill_with(short fill_value);
      void zeros(void) { fill_with(0); }
      void ones(void)  { fill_with(1); }

      MRI_Short_Matrix& operator=(const MRI_Short_Matrix& mat);
      MRI_Short_Matrix& operator*=(double a); 
      MRI_Short_Matrix& operator*=(const MRI_Short_Matrix& x);
      MRI_Short_Matrix& operator+=(const MRI_Short_Matrix& x);
      MRI_Short_Matrix& saxpy(double a, const MRI_Short_Matrix& x,
                              const MRI_Short_Matrix& y);
      MRI_Short_Matrix& gaxpy(const MRI_Short_Matrix& a, 
                              const MRI_Short_Matrix& x,
                              const MRI_Short_Matrix& y);

      int operator != (const MRI_Short_Matrix& mat) const;
      int operator == (const MRI_Short_Matrix& mat) const {
         return !(operator != (mat)); }

      short maximum(void) const;
      short minimum(void) const;
      double sum(void) const;
      double norm(void) const;
      double mean(void) const;
      double std(void) const;

      short maximum(unsigned int row1, unsigned int row2,
                             unsigned int col1, unsigned int col2) const;
      short minimum(unsigned int row1, unsigned int row2,
                             unsigned int col1, unsigned int col2) const;
      double sum(unsigned int row1, unsigned int row2,
                 unsigned int col1, unsigned int col2) const;
      double norm(unsigned int row1, unsigned int row2,
                  unsigned int col1, unsigned int col2) const;
      double mean(unsigned int row1, unsigned int row2,
                  unsigned int col1, unsigned int col2) const;
      double std(unsigned int row1, unsigned int row2,
                 unsigned int col1, unsigned int col2) const;

      void set_submatrix(MRI_Short_Matrix& mat, 
                         unsigned int row, 
                         unsigned int col);

      short *row_ptr(int row){
         return &(_matrix[row*_ncols]); }
      short *col_ptr(int col){
         return &(_matrix[col]); }

      const short *row_ptr(int row) const {
         return &(_matrix[row*_ncols]); }
      const short *col_ptr(int col) const {
         return &(_matrix[col]); }
   
      virtual operator void *() {
         return (void *)_matrix; }
      virtual operator const void *() const {
         return (const void *)_matrix; }

      operator short *() {
         return _matrix; }
      operator const short *() const {
         return _matrix; }

      short *operator[] (int row){
         return row_ptr(row); }
      const short *operator[] (int row) const {
         return row_ptr(row); }

      short& operator() (unsigned int row, unsigned int col){
         return _matrix[row*_ncols+col]; }
      short  operator() (unsigned int row, unsigned int col) const {
         return _matrix[row*_ncols+col]; }
      MRI_Short_Matrix operator() (unsigned int row1, unsigned int row2,
                                 unsigned int col1, unsigned int col2) const;

      void display(ostream& stream) const;

   protected:
      MRI_Short_Matrix();

      short *_matrix;

      const short *_matrix_endptr(void) const {
         return _matrix+get_nelements(); }

      virtual void _allocate(unsigned int nelements) {
         _matrix = new short[nelements]; }
      virtual void _deallocate(void) {
         delete[] _matrix; }

   private:     
      friend class MRI_Byte_Matrix;
      friend class MRI_Float_Matrix;
      friend class MRI_Double_Matrix;
      friend class MRI_FComplex_Matrix;
      friend class MRI_Complex_Matrix;
      friend class MRI_Label;
      friend class MRI_Image;
};

//===========================================================================
// MRI_Float_Matrix
//===========================================================================

class MRI_Float_Matrix : public MRI_Matrix {
   public:
      MRI_Float_Matrix(unsigned int nrows, unsigned int ncols, 
                       float fill = 0.0);
      MRI_Float_Matrix(const MRI_Float_Matrix& mat);
      MRI_Float_Matrix(const MRI_Byte_Matrix& mat);
      MRI_Float_Matrix(const MRI_Short_Matrix& mat);
      MRI_Float_Matrix(const MRI_Double_Matrix& mat);

      virtual ~MRI_Float_Matrix();

      virtual size_t size_in_bytes(void) const {
         return get_nelements()*sizeof(float); }
      virtual size_t element_size_in_bytes(void) const {
         return sizeof(float); }

      void fill_with(float fill_value);
      void zeros(void) { fill_with(0); }
      void ones(void)  { fill_with(1); }

      MRI_Float_Matrix& operator=(const MRI_Float_Matrix& mat);
      MRI_Float_Matrix& operator*=(double a);
      MRI_Float_Matrix& operator*=(const MRI_Float_Matrix& x);
      MRI_Float_Matrix& operator/=(const MRI_Float_Matrix& x);
      MRI_Float_Matrix& operator+=(const MRI_Float_Matrix& x);
      MRI_Float_Matrix& saxpy(double a, const MRI_Float_Matrix& x, 
                              const MRI_Float_Matrix& y);
      MRI_Float_Matrix& gaxpy(const MRI_Float_Matrix& a, 
                              const MRI_Float_Matrix& x, 
                              const MRI_Float_Matrix& y);

      int operator != (const MRI_Float_Matrix& mat) const;
      int operator == (const MRI_Float_Matrix& mat) const {
         return !(operator != (mat)); }

      double maximum(void) const;
      double minimum(void) const;
      double sum(void) const;
      double norm(void) const;
      double mean(void) const;
      double std(void) const;

      double maximum(unsigned int row1, unsigned int row2,
                     unsigned int col1, unsigned int col2) const;
      double minimum(unsigned int row1, unsigned int row2,
                     unsigned int col1, unsigned int col2) const;
      double sum(unsigned int row1, unsigned int row2,
                 unsigned int col1, unsigned int col2) const;
      double norm(unsigned int row1, unsigned int row2,
                  unsigned int col1, unsigned int col2) const;
      double mean(unsigned int row1, unsigned int row2,
                  unsigned int col1, unsigned int col2) const;
      double std(unsigned int row1, unsigned int row2,
                 unsigned int col1, unsigned int col2) const;

      MRI_FComplex_Matrix FFT(void) const;
      MRI_FComplex_Matrix iFFT(void) const;
      MRI_FComplex_Matrix FFT2(void) const;
      MRI_FComplex_Matrix iFFT2(void) const;
      void fftshift(void);

      void set_submatrix(MRI_Float_Matrix& mat, 
                         unsigned int row, unsigned int col);

      float *row_ptr(int row) {
         return &(_matrix[row*_ncols]); }
      float *col_ptr(int col) {
         return &(_matrix[col]); }

      const float *row_ptr(int row) const {
         return &(_matrix[row*_ncols]); }
      const float *col_ptr(int col) const {
         return &(_matrix[col]); }

      virtual operator void *() {
         return (void *)_matrix; }
      virtual operator const void *() const {
         return (const void *)_matrix; }

      operator float *() {
         return _matrix; }
      operator const float *() {
         return _matrix; }

      float *operator [] (int row) {
         return row_ptr(row); }
      const float *operator [] (int row) const {
         return row_ptr(row); }

      float& operator() (unsigned int row, unsigned int col){
         return _matrix[row*_ncols+col]; }
      float  operator() (unsigned int row, unsigned int col) const {
         return _matrix[row*_ncols+col]; }
      MRI_Float_Matrix operator() (unsigned int row1, unsigned int row2,
                                  unsigned int col1, unsigned int col2) const;

      void display(ostream& stream) const;

   protected:
      MRI_Float_Matrix();

      float *_matrix;
 
      const float *_matrix_endptr(void) const {
         return _matrix+get_nelements(); }
 
      virtual void _allocate(unsigned int nelements) {
         _matrix = new float[nelements]; }
      virtual void _deallocate(void) {
         delete[] _matrix; }

   private:
      friend class MRI_Byte_Matrix;
      friend class MRI_Short_Matrix;
      friend class MRI_Double_Matrix;
      friend class MRI_FComplex_Matrix;
      friend class MRI_Complex_Matrix;
      friend class MRI_Label;
      friend class MRI_Image;
};

//===========================================================================
// MRI_Double_Matrix
//===========================================================================

class MRI_Double_Matrix : public MRI_Matrix {
   public:
      MRI_Double_Matrix(unsigned int nrows, unsigned int ncols, 
                        double fill = 0.0);
      MRI_Double_Matrix(const MRI_Double_Matrix& mat);
      MRI_Double_Matrix(const MRI_Byte_Matrix& mat);
      MRI_Double_Matrix(const MRI_Short_Matrix& mat);
      MRI_Double_Matrix(const MRI_Float_Matrix& mat);

      virtual ~MRI_Double_Matrix();

      virtual size_t size_in_bytes(void) const {
         return get_nelements()*sizeof(double); }
      virtual size_t element_size_in_bytes(void) const {
         return sizeof(double); }

      void fill_with(double fill_value);
      void zeros(void) { fill_with(0); }
      void ones(void)  { fill_with(1); }

      MRI_Double_Matrix& operator=(const MRI_Double_Matrix& mat);
      MRI_Double_Matrix& operator*=(double a);
      MRI_Double_Matrix& operator*=(const MRI_Double_Matrix& x);
      MRI_Double_Matrix& operator/=(const MRI_Double_Matrix& x);
      MRI_Double_Matrix& operator+=(const MRI_Double_Matrix& x);
      MRI_Double_Matrix& saxpy(double a, const MRI_Double_Matrix& x,
                               const MRI_Double_Matrix& y);
      MRI_Double_Matrix& gaxpy(const MRI_Double_Matrix& a, 
                               const MRI_Double_Matrix& x,
                               const MRI_Double_Matrix& y);

      int operator != (const MRI_Double_Matrix& mat) const;
      int operator == (const MRI_Double_Matrix& mat) const {
         return !(operator != (mat)); }

      double maximum(void) const;
      double minimum(void) const;
      double sum(void) const;
      double norm(void) const;
      double mean(void) const;
      double std(void) const;

      double maximum(unsigned int row1, unsigned int row2,
                     unsigned int col1, unsigned int col2) const;
      double minimum(unsigned int row1, unsigned int row2,
                     unsigned int col1, unsigned int col2) const;
      double sum(unsigned int row1, unsigned int row2,
                 unsigned int col1, unsigned int col2) const;
      double norm(unsigned int row1, unsigned int row2,
                  unsigned int col1, unsigned int col2) const;
      double mean(unsigned int row1, unsigned int row2,
                  unsigned int col1, unsigned int col2) const;
      double std(unsigned int row1, unsigned int row2,
                 unsigned int col1, unsigned int col2) const;

      MRI_Complex_Matrix FFT(void) const;
      MRI_Complex_Matrix iFFT(void) const;
      MRI_Complex_Matrix FFT2(void) const;
      MRI_Complex_Matrix iFFT2(void) const;
      void fftshift(void);

      void set_submatrix(MRI_Double_Matrix& mat, 
                         unsigned int row, unsigned int col);

      double *row_ptr(int row) {
         return &(_matrix[row*_ncols]); }
      double *col_ptr(int col) {
         return &(_matrix[col]); }

      const double *row_ptr(int row) const {
         return &(_matrix[row*_ncols]); }
      const double *col_ptr(int col) const {
         return &(_matrix[col]); }

      virtual operator void *() {
         return (void *)_matrix; }
      virtual operator const void *() const {
         return (const void *)_matrix; }

      operator double *() {
         return _matrix; }
      operator const double *() const {
         return _matrix; }

      double *operator [] (int row){
         return row_ptr(row); }
      const double *operator [] (int row) const {
         return row_ptr(row); }

      double& operator() (unsigned int row, unsigned int col){
         return _matrix[row*_ncols+col]; }
      double  operator() (unsigned int row, unsigned int col) const {
         return _matrix[row*_ncols+col]; }
      MRI_Double_Matrix operator() (unsigned int row1, unsigned int row2,
                                  unsigned int col1, unsigned int col2) const;

      void display(ostream& stream) const;

   protected:
      MRI_Double_Matrix();

      double *_matrix;
 
      const double *_matrix_endptr(void) const {
         return _matrix+get_nelements(); }
 
      virtual void _allocate(unsigned int nelements) {
         _matrix = new double[nelements]; }
      virtual void _deallocate(void) {
         delete[] _matrix; }

   private:
      friend class MRI_Byte_Matrix;
      friend class MRI_Short_Matrix;
      friend class MRI_Float_Matrix;
      friend class MRI_FComplex_Matrix;
      friend class MRI_Complex_Matrix;
      friend class MRI_Label;
      friend class MRI_Image;
};

//===========================================================================
// MRI_FComplex_Matrix
//===========================================================================

class MRI_FComplex_Matrix : public MRI_Matrix {
   public:
      MRI_FComplex_Matrix(unsigned int nrows, 
                          unsigned int ncols, 
                          float fill = 0.0);
      MRI_FComplex_Matrix(const MRI_FComplex_Matrix& mat);
      MRI_FComplex_Matrix(const MRI_Byte_Matrix& mat);
      MRI_FComplex_Matrix(const MRI_Float_Matrix& mat);
      MRI_FComplex_Matrix(const MRI_Double_Matrix& mat);

      virtual ~MRI_FComplex_Matrix();

      virtual size_t size_in_bytes(void) const {
         return 2*get_nelements()*sizeof(float); }
      virtual size_t element_size_in_bytes(void) const {
         return 2*sizeof(float); }
      
      void fill_with(float fill_value);
      void fill_real_with(float fill_value);
      void fill_imag_with(float fill_value);
      void zeros(void) { fill_with(0); }
      void ones(void)  { fill_real_with(1); fill_imag_with(0);}

      void real(MRI_Float_Matrix& mat) const;
      void imag(MRI_Float_Matrix& mat) const;
      void abs(MRI_Float_Matrix& mat) const;
      void angle(MRI_Float_Matrix& mat) const;

      void get_real_min_max(double &min, double &max) const;
      void get_imag_min_max(double &min, double &max) const;
      void get_abs_min_max(double &min, double &max) const;
      void get_angle_min_max(double &min, double &max) const;

      MRI_FComplex_Matrix& operator=(const MRI_FComplex_Matrix& mat);
      MRI_FComplex_Matrix& operator*=(double a);
      MRI_FComplex_Matrix& operator*=(const MRI_FComplex_Matrix& x);
      MRI_FComplex_Matrix& operator*=(const MRI_Float_Matrix& x);
      MRI_FComplex_Matrix& operator/=(const MRI_Float_Matrix& x);
      MRI_FComplex_Matrix& operator+=(const MRI_FComplex_Matrix& x);
      MRI_FComplex_Matrix& saxpy(double areal, double aimag, 
                              const MRI_FComplex_Matrix& x,
                              const MRI_FComplex_Matrix& y);
      MRI_FComplex_Matrix& gaxpy(const MRI_FComplex_Matrix& a, 
                              const MRI_FComplex_Matrix& x,
                              const MRI_FComplex_Matrix& y);

      void FFT(void);
      void iFFT(void);
      void FFT2(void);
      void iFFT2(void);
      void fftshift(void);

      void set_submatrix(MRI_FComplex_Matrix& mat, 
                         unsigned int row, unsigned int col);
      void set_submatrix(MRI_Float_Matrix& mat, 
                         unsigned int row, unsigned int col);

      float *row_ptr(int row) {
         return &(_matrix[2*row*_ncols]); }
      float *col_ptr(int col) {
         return &(_matrix[2*col]); }

      const float *row_ptr(int row) const {
         return &(_matrix[2*row*_ncols]); }
      const float *col_ptr(int col) const {
         return &(_matrix[2*col]); }

      virtual operator void *() {
         return (void *)_matrix; }
      virtual operator const void *() const {
         return (const void *)_matrix; }

      operator float *() {
         return _matrix; }
      operator const float *() const {
         return _matrix; }

      float *operator[](int row) {
         return row_ptr(row); }
      const float *operator[](int row) const {
         return row_ptr(row); }

      float& real(unsigned int row, unsigned int col){
         return _matrix[2*(row*_ncols+col)]; }
      float& imag(unsigned int row, unsigned int col){
         return _matrix[2*(row*_ncols+col)+1]; }

      float real(unsigned int row, unsigned int col) const {
         return _matrix[2*(row*_ncols+col)]; }
      float imag(unsigned int row, unsigned int col) const {
         return _matrix[2*(row*_ncols+col)+1]; }
      float  abs(unsigned int row, unsigned int col) const {
         return hypotf(real(row,col), imag(row,col)); }
      float  angle(unsigned int row, unsigned int col) const {
         return (float) atan2(real(row,col), imag(row,col)); }

      void display(ostream& stream) const;

   protected:
      MRI_FComplex_Matrix();

      float *_matrix;

      const float *_matrix_endptr(void) const {
         return _matrix+2*get_nelements(); }

      virtual void _allocate(unsigned int nelements) {
         _matrix = new float[2*nelements]; }
      virtual void _deallocate(void) {
         delete[] _matrix; }

   private:
      friend class MRI_Complex_Matrix;

};

//===========================================================================
// MRI_Complex_Matrix
//===========================================================================

class MRI_Complex_Matrix : public MRI_Matrix {
   public:
      MRI_Complex_Matrix(unsigned int nrows, 
                         unsigned int ncols, 
                         double fill = 0.0);
      MRI_Complex_Matrix(const MRI_Complex_Matrix& mat);
      MRI_Complex_Matrix(const MRI_FComplex_Matrix& mat);      
      MRI_Complex_Matrix(const MRI_Byte_Matrix& mat);      
      MRI_Complex_Matrix(const MRI_Float_Matrix& mat);
      MRI_Complex_Matrix(const MRI_Double_Matrix& mat);

      virtual ~MRI_Complex_Matrix();

      virtual size_t size_in_bytes(void) const {
         return 2*get_nelements()*sizeof(double); }
      virtual size_t element_size_in_bytes(void) const {
         return 2*sizeof(double); }
      
      void fill_with(double fill_value);
      void fill_real_with(double fill_value);
      void fill_imag_with(double fill_value);
      void zeros(void) { fill_with(0); }
      void ones(void)  { fill_real_with(1); fill_imag_with(0); }

      void real(MRI_Double_Matrix& mat) const;
      void imag(MRI_Double_Matrix& mat) const;
      void abs(MRI_Double_Matrix& mat) const;
      void angle(MRI_Double_Matrix& mat) const;

      void get_real_min_max(double &min, double &max) const;
      void get_imag_min_max(double &min, double &max) const;
      void get_abs_min_max(double &min, double &max) const;
      void get_angle_min_max(double &min, double &max) const;

      MRI_Complex_Matrix& operator=(const MRI_Complex_Matrix& mat);
      MRI_Complex_Matrix& operator*=(double a);
      MRI_Complex_Matrix& operator*=(const MRI_Complex_Matrix& x);
      MRI_Complex_Matrix& operator*=(const MRI_Double_Matrix& x);
      MRI_Complex_Matrix& operator/=(const MRI_Double_Matrix& x);
      MRI_Complex_Matrix& operator+=(const MRI_Complex_Matrix& x);
      MRI_Complex_Matrix& saxpy(double areal, double aimag, 
                              const MRI_Complex_Matrix& x,
                              const MRI_Complex_Matrix& y);
      MRI_Complex_Matrix& gaxpy(const MRI_Complex_Matrix& a, 
                              const MRI_Complex_Matrix& x,
                              const MRI_Complex_Matrix& y);

      void FFT(void);
      void iFFT(void);
      void FFT2(void);
      void iFFT2(void);
      void fftshift(void);

      void set_submatrix(MRI_Complex_Matrix& mat, 
                         unsigned int row, unsigned int col);
      void set_submatrix(MRI_Byte_Matrix& mat, 
                         unsigned int row, unsigned int col);

      double *row_ptr(int row) {
         return &(_matrix[2*row*_ncols]); }
      double *col_ptr(int col) {
         return &(_matrix[2*col]); }

      const double *row_ptr(int row) const {
         return &(_matrix[2*row*_ncols]); }
      const double *col_ptr(int col) const {
         return &(_matrix[2*col]); }

      virtual operator void *() {
         return (void *)_matrix; }
      virtual operator const void *() const {
         return (const void *)_matrix; }

      operator double *() {
         return _matrix; }
      operator const double *() const {
         return _matrix; }

      double *operator[](int row) {
         return row_ptr(row); }
      const double *operator[](int row) const {
         return row_ptr(row); }

      double& real(unsigned int row, unsigned int col){
         return _matrix[2*(row*_ncols+col)]; }
      double& imag(unsigned int row, unsigned int col){
         return _matrix[2*(row*_ncols+col)+1]; }

      double  real(unsigned int row, unsigned int col) const {
         return _matrix[2*(row*_ncols+col)]; }
      double  imag(unsigned int row, unsigned int col) const {
         return _matrix[2*(row*_ncols+col)+1]; }
      double  abs(unsigned int row, unsigned int col) const {
         return hypot(real(row,col), imag(row,col)); }
      double  angle(unsigned int row, unsigned int col) const {
         return atan2(real(row,col), imag(row,col)); }

      void display(ostream& stream) const;

   protected:
      MRI_Complex_Matrix();

      double *_matrix;

      const double *_matrix_endptr(void) const {
         return _matrix+2*get_nelements(); }

      virtual void _allocate(unsigned int nelements) {
         _matrix = new double[2*nelements]; }
      virtual void _deallocate(void) {
         delete[] _matrix; }


};

#endif





