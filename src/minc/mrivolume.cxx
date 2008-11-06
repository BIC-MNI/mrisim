//===========================================================================
// MRIVOLUME.CXX
//
// R.Kwan
// August 31, 1995
//
// (C) Copyright 1995 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/mrivolume.cxx,v 1.2 2008-11-06 10:58:23 rotor Exp $
// $Log: mrivolume.cxx,v $
// Revision 1.2  2008-11-06 10:58:23  rotor
//  * fixed includes for iostream and friends
//  * updated for new release (1.0.2)
//
// Revision 1.1  2003/05/30 16:43:09  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.6  1996/05/29  19:05:28  rkwan
// GNU warning fix
//
// Revision 2.5  1996/05/29  16:29:01  rkwan
// Release 2.5
//
// Revision 1.6  1996/01/08  15:34:16  rkwan
// Added MRI_Float_Volume.
//
// Revision 1.5  1995/12/11  16:13:04  rkwan
// fix MRI_Short_Volume::display
//
// Revision 1.4  1995/12/11  14:18:12  rkwan
// Updated for fast_iso_model.
//
//===========================================================================

#include "mrivolume.h"

#include <math.h>
#include <iostream>
#include <iomanip>
#include <limits.h>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

//===========================================================================
// MRI_Volume
//===========================================================================

//---------------------------------------------------------------------------
// MRI_Volume constructors
//---------------------------------------------------------------------------

MRI_Volume::MRI_Volume() {
   _nrows = 0;
   _ncols = 0;
   _nslices = 0;
}

MRI_Volume::MRI_Volume(unsigned int nrows, unsigned int ncols,
                     unsigned int nslices) {
   _nrows = nrows;
   _ncols = ncols;
   _nslices = nslices;
}

//---------------------------------------------------------------------------
// MRI_Volume destructor
//---------------------------------------------------------------------------

MRI_Volume::~MRI_Volume() {}

//---------------------------------------------------------------------------
// MRI_Volume::is_same_size_as
// Returns TRUE if the two matrices have the same dimensions.
//---------------------------------------------------------------------------

int MRI_Volume::is_same_size_as(MRI_Volume& vol) const {
   return (((_nrows == vol._nrows) &&
            (_ncols == vol._ncols) &&
            (_nslices == vol._nslices)) ? TRUE : FALSE);
}

//---------------------------------------------------------------------------
// MRI_Volume::is_same_size_as
// Returns TRUE if the two matrices have the same dimensions.
//---------------------------------------------------------------------------

int MRI_Volume::slice_size_is_same_as(MRI_Matrix& mat) const {
   return (((_nrows == mat.get_nrows()) &&
            (_ncols == mat.get_ncols())) ? TRUE : FALSE);
}

//---------------------------------------------------------------------------
// MRI_Volume::is_power_of_two
// Returns TRUE if the matrix dimensions are a power of two.
//---------------------------------------------------------------------------

int MRI_Volume::is_power_of_two(void) const {
   return ((_power_of_two(_nrows) &&
            _power_of_two(_ncols) &&
            _power_of_two(_nslices)) ? TRUE : FALSE);
}

//============================================================================
// MRI_Byte_Volume
//============================================================================

//---------------------------------------------------------------------------
// MRI_Byte_Volume constructors
//---------------------------------------------------------------------------

MRI_Byte_Volume::MRI_Byte_Volume() : MRI_Volume() {
   _allocate(this->get_nelements());
}

MRI_Byte_Volume::MRI_Byte_Volume(unsigned int nrows, 
                             unsigned int ncols,
                             unsigned int nslices,
                             unsigned char fill_value) 
   : MRI_Volume(nrows, ncols, nslices) {

   _allocate(this->get_nelements());
   this->fill_with(fill_value);
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume copy constructor 
//---------------------------------------------------------------------------

MRI_Byte_Volume::MRI_Byte_Volume(const MRI_Byte_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   _allocate(this->get_nelements());
   (void)memcpy(_volume, vol._volume, this->size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume type cast constructor
//---------------------------------------------------------------------------

MRI_Byte_Volume::MRI_Byte_Volume(const MRI_Short_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);
   for(n=0; n<len; n++){
      _volume[n] = (short) (vol._volume[n]);
   }
}

MRI_Byte_Volume::MRI_Byte_Volume(const MRI_Float_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);
   for(n=0; n<len; n++){
      _volume[n] = (unsigned char) rint(vol._volume[n]);
   }
}

MRI_Byte_Volume::MRI_Byte_Volume(const MRI_Double_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);
   for(n=0; n<len; n++){
      _volume[n] = (unsigned char) rint(vol._volume[n]);
   }
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume destructor
//---------------------------------------------------------------------------

MRI_Byte_Volume::~MRI_Byte_Volume() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Byte_Volume::fill_with(unsigned char fill_value){
   (void)memset(_volume, fill_value, this->get_nelements());
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume:operator=
// Assignment operator
//---------------------------------------------------------------------------

MRI_Byte_Volume& MRI_Byte_Volume::operator=(const MRI_Byte_Volume& vol){
   if (this != &vol){
      _nrows = vol._nrows;
      _ncols = vol._ncols;
      _nslices = vol._nslices;
      (void)memcpy(_volume, vol._volume, this->get_nelements());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume::operator*=
// Scales a matrix by a (double) scalar.
//---------------------------------------------------------------------------

MRI_Byte_Volume& MRI_Byte_Volume::operator*=(double a) { 
   unsigned char *ptr;
   for(ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      *ptr = (unsigned char)rint(*ptr * a);
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume::operator !=
// Element-wise comparison.
//---------------------------------------------------------------------------

int MRI_Byte_Volume::operator != (const MRI_Byte_Volume& vol) const {
   unsigned char *source, *target;
   int equal;
   for(equal=1, source=vol._volume, target=_volume;
       (equal==1) && target<_volume+get_nelements(); 
       source++, target++){
     if (*source != *target) equal = 0;
   }
   return equal;
}  

//---------------------------------------------------------------------------
// MRI_Byte_Volume::maximum 
// Returns the maximum element in the matrix.
//---------------------------------------------------------------------------

unsigned char MRI_Byte_Volume::maximum(void) const {
   unsigned char max, *ptr;
   for(max=0, ptr=_volume;
       ptr<_volume+get_nelements();
       ptr++){
      if (*ptr > max) max = *ptr;
   }
   return max;
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume::minimum
// Returns the minimum element in the matrix.
//---------------------------------------------------------------------------

unsigned char MRI_Byte_Volume::minimum(void) const{
   unsigned char min, *ptr;
   for(min=CHAR_MAX, ptr=_volume;
       ptr<_volume+get_nelements();
       ptr++){
      if (*ptr < min) min = *ptr;
   }
   return min;
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume::sum
// Returns the sum of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Byte_Volume::sum(void) const {
   double sum;
   unsigned char *ptr;
   for(sum=0.0,ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      sum += (double)*ptr;
   }
   return sum;
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume::norm
// Returns the norm of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Byte_Volume::norm(void) const{
   double avg = this->mean();
   double sum;
   unsigned char *ptr;

   for(sum=0.0,ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      sum += SQR((double)*ptr-avg);
   }
   return sqrt(sum);
}
 
//---------------------------------------------------------------------------
// MRI_Byte_Volume::mean
// Returns the mean of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Byte_Volume::mean(void) const{
   double sum = this->sum(); 
   return sum/get_nelements();
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume::std
// Returns the standard deviation of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Byte_Volume::std(void) const{
   double norm = this->norm();
   return norm/sqrt(get_nelements()-1);
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume::operator () 
// Returns a submatrix of the matrix.
//---------------------------------------------------------------------------

MRI_Byte_Volume MRI_Byte_Volume::operator() (unsigned int row1, 
                                         unsigned int row2,
                                         unsigned int col1, 
                                         unsigned int col2,
                                         unsigned int slice1, 
                                         unsigned int slice2) const {

#ifdef DEBUG
   assert(row2 >= row1); 
   assert(col2 >= col1);
   assert(slice2 >= slice1);
#endif

   unsigned int row, col, slice;
   MRI_Byte_Volume vol(row2-row1+1,col2-col1+1,slice2-slice1+1);
   for(slice=slice1; slice<=slice2; slice++){
      for(row=row1; row<=row2; row++){
         for(col=col1; col<=col2; col++){
            vol(row-row1,col-col1,slice-slice1) = (*this)(row,col,slice);
         }
      }
   }
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Byte_Volume::display
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_Byte_Volume::display(ostream& stream){
   unsigned int row, col, slice;
   for(slice=0; slice<_nslices; slice++){
      stream << endl;
      for(row=0; row<_nrows; row++){
         for(col=0; col<_ncols; col++){
            stream << (int)(*this)(row,col,slice) << " ";
         }
         stream << endl;
      }
   }
}

//============================================================================
// MRI_Short_Volume
//============================================================================

//---------------------------------------------------------------------------
// MRI_Short_Volume constructors
//---------------------------------------------------------------------------

MRI_Short_Volume::MRI_Short_Volume() : MRI_Volume() {
   _allocate(this->get_nelements());
}

MRI_Short_Volume::MRI_Short_Volume(unsigned int nrows,
                               unsigned int ncols,
                               unsigned int nslices,
                               short fill_value)
   : MRI_Volume(nrows, ncols, nslices) {

   _allocate(this->get_nelements());
   this->fill_with(fill_value);
}

//---------------------------------------------------------------------------
// MRI_Short_Volume copy constructor
//---------------------------------------------------------------------------

MRI_Short_Volume::MRI_Short_Volume(const MRI_Short_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   _allocate(this->get_nelements());
   (void)memcpy(_volume, vol._volume, this->size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Short_Volume type cast constructor
//---------------------------------------------------------------------------

MRI_Short_Volume::MRI_Short_Volume(const MRI_Byte_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   unsigned char  *source;
   short *target;

   _allocate(this->get_nelements());
   for (source=vol._volume, target=_volume;
        target<_volume+get_nelements();
        source++, target++){

      *target = (short) rint(*source);
   }
}

MRI_Short_Volume::MRI_Short_Volume(const MRI_Float_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   float *source;
   short *target;

   _allocate(this->get_nelements());
   for (source=vol._volume, target=_volume;
        target<_volume+get_nelements();
        source++, target++){

      *target = (short) rint(*source);
   }
}

MRI_Short_Volume::MRI_Short_Volume(const MRI_Double_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   double        *source;
   short *target;

   _allocate(this->get_nelements());
   for (source=vol._volume, target=_volume;
        target<_volume+get_nelements();
        source++, target++){

      *target = (short) rint(*source);
   }
}

//---------------------------------------------------------------------------
// MRI_Short_Volume destructor
//---------------------------------------------------------------------------

MRI_Short_Volume::~MRI_Short_Volume() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Short_Volume::fill_with(short fill_value){
   (void)memset(_volume, fill_value, this->get_nelements());
}

//---------------------------------------------------------------------------
// MRI_Short_Volume:operator=
// Assignment operator
//---------------------------------------------------------------------------

MRI_Short_Volume& MRI_Short_Volume::operator=(const MRI_Short_Volume& vol){
   if (this != &vol){
      _nrows = vol._nrows;
      _ncols = vol._ncols;
      _nslices = vol._nslices;
      (void)memcpy(_volume, vol._volume, this->get_nelements());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::operator*=
// Scales a matrix by a (double) scalar.
//---------------------------------------------------------------------------

MRI_Short_Volume& MRI_Short_Volume::operator*=(double a) {
   short *ptr;
   for(ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      *ptr = (short)rint(*ptr * a);
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::operator !=
// Element-wise comparison.
//---------------------------------------------------------------------------

int MRI_Short_Volume::operator != (const MRI_Short_Volume& vol) const {
   short *source, *target;
   int equal;
   for(equal=1, source=vol._volume, target=_volume;
       (equal==1) && target<_volume+get_nelements();
       source++, target++){
     if (*source != *target) equal = 0;
   }
   return equal;
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::maximum
// Returns the maximum element in the matrix.
//---------------------------------------------------------------------------

short MRI_Short_Volume::maximum(void) const {
   short max, *ptr;
   for(max=SHRT_MIN, ptr=_volume;
       ptr<_volume+get_nelements();
       ptr++){
      if (*ptr > max) max = *ptr;
   }
   return max;
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::minimum
// Returns the minimum element in the matrix.
//---------------------------------------------------------------------------

short MRI_Short_Volume::minimum(void) const{
   short min, *ptr;
   for(min=SHRT_MAX, ptr=_volume;
       ptr<_volume+get_nelements();
       ptr++){
      if (*ptr < min) min = *ptr;
   }
   return min;
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::sum
// Returns the sum of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Short_Volume::sum(void) const {
   double sum;
   short *ptr;
   for(sum=0.0,ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      sum += (double)*ptr;
   }
   return sum;
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::norm
// Returns the norm of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Short_Volume::norm(void) const{
   double avg = this->mean();
   double sum;
   short *ptr;

   for(sum=0.0,ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      sum += SQR((double)*ptr-avg);
   }
   return sqrt(sum);
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::mean
// Returns the mean of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Short_Volume::mean(void) const{
   double sum = this->sum();
   return sum/get_nelements();
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::std
// Returns the standard deviation of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Short_Volume::std(void) const{
   double norm = this->norm();
   return norm/sqrt(get_nelements()-1);
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::operator ()
// Returns a submatrix of the matrix.
//---------------------------------------------------------------------------

MRI_Short_Volume MRI_Short_Volume::operator() (unsigned int row1, 
                                           unsigned int row2,
                                           unsigned int col1, 
                                           unsigned int col2,
                                           unsigned int slice1, 
                                           unsigned int slice2) const {

#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
   assert(slice2 >= slice1);
#endif

   unsigned int row, col, slice;
   MRI_Short_Volume vol(row2-row1+1,col2-col1+1,slice2-slice1+1);
   for(slice=slice1; slice<=slice2; slice++){
      for(row=row1; row<=row2; row++){
         for(col=col1; col<=col2; col++){
            vol(row-row1,col-col1,slice-slice1) = (*this)(row,col,slice);
         }
      }
   }
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Short_Volume::display
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_Short_Volume::display(ostream& stream){
   unsigned int row, col, slice;
   for(slice=0; slice<_nslices; slice++){
      stream << endl;
      for(row=0; row<_nrows; row++){
         for(col=0; col<_ncols; col++){
            stream << (short)(*this)(row,col,slice) << " ";
         }
         stream << endl;
      }
   }
}

//===========================================================================
// MRI_Float_Volume
//===========================================================================

//---------------------------------------------------------------------------
// MRI_Float_Volume constructors
//---------------------------------------------------------------------------

MRI_Float_Volume::MRI_Float_Volume() : MRI_Volume() {
   _allocate(this->get_nelements());
}

MRI_Float_Volume::MRI_Float_Volume(unsigned int nrows, 
                                 unsigned int ncols,
                                 unsigned int nslices,
                                 float fill_value) :
    MRI_Volume(nrows, ncols, nslices) {

   _allocate(this->get_nelements());
   this->fill_with(fill_value);
}

//---------------------------------------------------------------------------
// MRI_Float_Volume copy constructor
//---------------------------------------------------------------------------

MRI_Float_Volume::MRI_Float_Volume(const MRI_Float_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   _allocate(this->get_nelements());
   (void)memcpy(_volume, vol._volume, this->size_in_bytes());

}

//---------------------------------------------------------------------------
// MRI_Float_Volume type cast constructor
//---------------------------------------------------------------------------

MRI_Float_Volume::MRI_Float_Volume(const MRI_Byte_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   _allocate(this->get_nelements());
   float *target;
   unsigned char *source;

   for (source=vol._volume, target=_volume;
        target<_volume+get_nelements();
        source++, target++){
      *target = (float)*source;
   }
}
 
MRI_Float_Volume::MRI_Float_Volume(const MRI_Short_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   _allocate(this->get_nelements());
   float *target;
   short *source;

   for (source=vol._volume, target=_volume;
        target<_volume+get_nelements();
        source++, target++){
      *target = (float)*source;
   }
}

MRI_Float_Volume::MRI_Float_Volume(const MRI_Double_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   _allocate(this->get_nelements());
   float *target;
   double *source;

   for (source=vol._volume, target=_volume;
        target<_volume+get_nelements();
        source++, target++){
      *target = (float)*source;
   }
}

//---------------------------------------------------------------------------
// MRI_Float_Volume destructor
//---------------------------------------------------------------------------

MRI_Float_Volume::~MRI_Float_Volume() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Float_Volume::fill_with(float fill_value) {

   float *ptr; 
   for (ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      *ptr = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::operator=
// Assignment operator.
//---------------------------------------------------------------------------

MRI_Float_Volume& MRI_Float_Volume::operator=(const MRI_Float_Volume& vol){
   if (this != &vol){
      _nrows = vol._nrows;
      _ncols = vol._ncols;
      _nslices = vol._nslices;
      (void)memcpy(_volume, vol._volume, this->size_in_bytes());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::operator*=
// Scales a matrix by a (float) scalar.
//---------------------------------------------------------------------------

MRI_Float_Volume& MRI_Float_Volume::operator*=(double a){
   float *ptr;
   for(ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      *ptr *= a;
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::operator!=
// Element-wise comparison.
//---------------------------------------------------------------------------

int MRI_Float_Volume::operator != (const MRI_Float_Volume& vol) const {
   float *source, *target;
   int equal;
   for (equal=1, source=vol._volume, target=_volume;
        (equal==1) && target<_volume+get_nelements();
        source++, target++){
      if (*source != *target) equal = 0;
   }
   return equal;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::maximum
// Returns the maximum element in the matrix.
//---------------------------------------------------------------------------

float MRI_Float_Volume::maximum(void) const {
   float max, *ptr;
   for(max=FLT_MIN, ptr=_volume;
       ptr<_volume+get_nelements();
       ptr++){
      if  (*ptr > max) max = *ptr;
   }
   return max;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::minimum
// Returns the minimum element in the matrix.
//---------------------------------------------------------------------------

float MRI_Float_Volume::minimum(void) const {
   float min, *ptr;
   for(min=FLT_MAX, ptr=_volume;
       ptr<_volume+get_nelements();
       ptr++){
      if (*ptr < min) min = *ptr;
   }
   return min;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::sum
// Returns the sum of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Float_Volume::sum(void) const {
   double sum;
   float  *ptr;
   for(sum=0.0, ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      sum += *ptr;
   }
   return sum;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::norm
// Returns the norm of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Float_Volume::norm(void) const {
   double avg = this->mean();
   double sum;
   float *ptr;
   for(sum=0.0, ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      sum += SQR(*ptr-avg);
   }
   return sqrt(sum);
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::mean
// Returns the mean of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Float_Volume::mean(void) const {
   double sum = this->sum();
   return sum/get_nelements();
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::std
// Returns the standard deviation of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Float_Volume::std(void) const {
   double norm = this->norm();
   return norm/sqrt(get_nelements()-1);
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::FFT2
// Returns the 2D-FFT of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Volume MRI_Float_Volume::FFT2(void) const {
   MRI_Complex_Volume vol(*this); 
   vol.FFT2();
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::iFFT2
// Returns the 2D-inverse FFT of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Volume MRI_Float_Volume::iFFT2(void) const {
   MRI_Complex_Volume vol(*this);
   vol.iFFT2();
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::FFT3
// Returns the 3D-FFT of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Volume MRI_Float_Volume::FFT3(void) const {
   MRI_Complex_Volume vol(*this);
   vol.FFT3();
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::iFFT3
// Returns the 3D-inverse FFT of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Volume MRI_Float_Volume::iFFT3(void) const {
   MRI_Complex_Volume vol(*this);
   vol.iFFT3();
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::fftshift2
// Swaps first and fourth, second and third quadrants to move the
// zeroth lag to the centre of the spectrum.
//---------------------------------------------------------------------------

void MRI_Float_Volume::fftshift2(void){
   unsigned int islice;
   for(islice=0; islice<get_nslices(); islice++){
      _fftshift((void *)(*this)[islice], _nrows, _ncols, this->element_size_in_bytes());
   }
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::fftshift3
// Swaps first and fourth, second and third quadrants to move the
// zeroth lag to the centre of the spectrum.
//---------------------------------------------------------------------------

void MRI_Float_Volume::fftshift3(void){
// Swap octants
//   _fftshift((void *)_volume, _nrows, _ncols, this->element_size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::operator()
// Returns a submatrix of the matrix.
//---------------------------------------------------------------------------

MRI_Float_Volume MRI_Float_Volume::operator() (unsigned int row1, unsigned int row2,
                                 unsigned int col1, unsigned int col2,
                                 unsigned int slice1, unsigned int slice2) const {

#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
   assert(slice2 >= slice1);
#endif

   unsigned int row, col, slice;
   MRI_Float_Volume vol(row2-row1+1,col2-col1+1,slice2-slice1+1);
   for(slice=slice1; slice<=slice2; slice++){
      for(row=row1; row<=row2; row++){
         for(col=col1; col<=col2; col++){
            vol(row-row1,col-col1,slice-slice1) = (*this)(row,col,slice);
         }
      }
   }
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Float_Volume::display
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_Float_Volume::display(ostream& stream){
   unsigned int row, col, slice;
   for(slice=0; slice<_nslices; slice++){
      stream << endl;
      for(row=0; row<_nrows; row++){
         for(col=0; col<_ncols; col++){
            stream << setprecision(4) << (*this)(row,col,slice) << " ";
         }
         stream << endl;
      }
   }
}

//===========================================================================
// MRI_Double_Volume
//===========================================================================

//---------------------------------------------------------------------------
// MRI_Double_Volume constructors
//---------------------------------------------------------------------------

MRI_Double_Volume::MRI_Double_Volume() : MRI_Volume() {
   _allocate(this->get_nelements());
}

MRI_Double_Volume::MRI_Double_Volume(unsigned int nrows, 
                                 unsigned int ncols,
                                 unsigned int nslices,
                                 double fill_value) :
    MRI_Volume(nrows, ncols, nslices) {

   _allocate(this->get_nelements());
   this->fill_with(fill_value);
}

//---------------------------------------------------------------------------
// MRI_Double_Volume copy constructor
//---------------------------------------------------------------------------

MRI_Double_Volume::MRI_Double_Volume(const MRI_Double_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   _allocate(this->get_nelements());
   (void)memcpy(_volume, vol._volume, this->size_in_bytes());

}

//---------------------------------------------------------------------------
// MRI_Double_Volume type cast constructor
//---------------------------------------------------------------------------

MRI_Double_Volume::MRI_Double_Volume(const MRI_Byte_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   _allocate(this->get_nelements());
   double *target;
   unsigned char *source;

   for (source=vol._volume, target=_volume;
        target<_volume+get_nelements();
        source++, target++){
      *target = (double)*source;
   }
}
 
MRI_Double_Volume::MRI_Double_Volume(const MRI_Short_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   _allocate(this->get_nelements());
   double *target;
   short *source;

   for (source=vol._volume, target=_volume;
        target<_volume+get_nelements();
        source++, target++){
      *target = (double)*source;
   }
}

//---------------------------------------------------------------------------
// MRI_Double_Volume destructor
//---------------------------------------------------------------------------

MRI_Double_Volume::~MRI_Double_Volume() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Double_Volume::fill_with(double fill_value) {

   double *ptr; 
   for (ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      *ptr = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::operator=
// Assignment operator.
//---------------------------------------------------------------------------

MRI_Double_Volume& MRI_Double_Volume::operator=(const MRI_Double_Volume& vol){
   if (this != &vol){
      _nrows = vol._nrows;
      _ncols = vol._ncols;
      _nslices = vol._nslices;
      (void)memcpy(_volume, vol._volume, this->size_in_bytes());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::operator*=
// Scales a matrix by a (double) scalar.
//---------------------------------------------------------------------------

MRI_Double_Volume& MRI_Double_Volume::operator*=(double a){
   double *ptr;
   for(ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      *ptr *= a;
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::operator!=
// Element-wise comparison.
//---------------------------------------------------------------------------

int MRI_Double_Volume::operator != (const MRI_Double_Volume& vol) const {
   double *source, *target;
   int equal;
   for (equal=1, source=vol._volume, target=_volume;
        (equal==1) && target<_volume+get_nelements();
        source++, target++){
      if (*source != *target) equal = 0;
   }
   return equal;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::maximum
// Returns the maximum element in the matrix.
//---------------------------------------------------------------------------

double MRI_Double_Volume::maximum(void) const {
   double max, *ptr;
   for(max=DBL_MIN, ptr=_volume;
       ptr<_volume+get_nelements();
       ptr++){
      if (*ptr > max) max = *ptr;
   }
   return max;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::minimum
// Returns the minimum element in the matrix.
//---------------------------------------------------------------------------

double MRI_Double_Volume::minimum(void) const {
   double min, *ptr;
   for(min=DBL_MAX, ptr=_volume;
       ptr<_volume+get_nelements();
       ptr++){
      if (*ptr < min) min = *ptr;
   }
   return min;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::sum
// Returns the sum of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Double_Volume::sum(void) const {
   double sum, *ptr;
   for(sum=0.0, ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      sum += *ptr;
   }
   return sum;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::norm
// Returns the norm of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Double_Volume::norm(void) const {
   double avg = this->mean();
   double sum, *ptr;
   for(sum=0.0, ptr=_volume; ptr<_volume+get_nelements(); ptr++){
      sum += SQR(*ptr-avg);
   }
   return sqrt(sum);
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::mean
// Returns the mean of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Double_Volume::mean(void) const {
   double sum = this->sum();
   return sum/get_nelements();
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::std
// Returns the standard deviation of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Double_Volume::std(void) const {
   double norm = this->norm();
   return norm/sqrt(get_nelements()-1);
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::FFT2
// Returns the 2D-FFT of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Volume MRI_Double_Volume::FFT2(void) const {
   MRI_Complex_Volume vol(*this); 
   vol.FFT2();
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::iFFT2
// Returns the 2D-inverse FFT of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Volume MRI_Double_Volume::iFFT2(void) const {
   MRI_Complex_Volume vol(*this);
   vol.iFFT2();
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::FFT3
// Returns the 3D-FFT of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Volume MRI_Double_Volume::FFT3(void) const {
   MRI_Complex_Volume vol(*this);
   vol.FFT3();
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::iFFT3
// Returns the 3D-inverse FFT of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Volume MRI_Double_Volume::iFFT3(void) const {
   MRI_Complex_Volume vol(*this);
   vol.iFFT3();
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::fftshift2
// Swaps first and fourth, second and third quadrants to move the
// zeroth lag to the centre of the spectrum.
//---------------------------------------------------------------------------

void MRI_Double_Volume::fftshift2(void){
   unsigned int islice;
   for(islice=0; islice<get_nslices(); islice++){
      _fftshift((void *)(*this)[islice], _nrows, _ncols, this->element_size_in_bytes());
   }
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::fftshift3
// Swaps first and fourth, second and third quadrants to move the
// zeroth lag to the centre of the spectrum.
//---------------------------------------------------------------------------

void MRI_Double_Volume::fftshift3(void){
// Swap octants
//   _fftshift((void *)_volume, _nrows, _ncols, this->element_size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::operator()
// Returns a submatrix of the matrix.
//---------------------------------------------------------------------------

MRI_Double_Volume MRI_Double_Volume::operator() (unsigned int row1, unsigned int row2,
                                 unsigned int col1, unsigned int col2,
                                 unsigned int slice1, unsigned int slice2) const {

#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
   assert(slice2 >= slice1);
#endif

   unsigned int row, col, slice;
   MRI_Double_Volume vol(row2-row1+1,col2-col1+1,slice2-slice1+1);
   for(slice=slice1; slice<=slice2; slice++){
      for(row=row1; row<=row2; row++){
         for(col=col1; col<=col2; col++){
            vol(row-row1,col-col1,slice-slice1) = (*this)(row,col,slice);
         }
      }
   }
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::display
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_Double_Volume::display(ostream& stream){
   unsigned int row, col, slice;
   for(slice=0; slice<_nslices; slice++){
      stream << endl;
      for(row=0; row<_nrows; row++){
         for(col=0; col<_ncols; col++){
            stream << setprecision(4) << (*this)(row,col,slice) << " ";
         }
         stream << endl;
      }
   }
}

//===========================================================================
// MRI_Complex_Volume
//===========================================================================

//---------------------------------------------------------------------------
// MRI_Complex_Volume constructors
//---------------------------------------------------------------------------

MRI_Complex_Volume::MRI_Complex_Volume() : MRI_Volume() {
   _allocate(2*this->get_nelements());
}

MRI_Complex_Volume::MRI_Complex_Volume(unsigned int nrows, 
                                   unsigned int ncols,
                                   unsigned int nslices,
                                   double fill_value) :
    MRI_Volume(nrows, ncols, nslices) {

   _allocate(2*this->get_nelements());
   this->fill_with(fill_value);
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume copy constructor
//---------------------------------------------------------------------------

MRI_Complex_Volume::MRI_Complex_Volume(const MRI_Complex_Volume& vol) :
    MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   _allocate(2*this->get_nelements());
   (void)memcpy(_volume, vol._volume, this->size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume type cast constructor
//---------------------------------------------------------------------------

MRI_Complex_Volume::MRI_Complex_Volume(const MRI_Float_Volume& vol) :
   MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   float *source;
   double *target;
   _allocate(2*this->get_nelements());
   for (source=vol._volume, target=_volume;
        target < _volume+2*get_nelements();
        source++, target += 2){
      *target     = (double)*source;
      *(target+1) = 0.0;
   }
}

MRI_Complex_Volume::MRI_Complex_Volume(const MRI_Double_Volume& vol) :
   MRI_Volume(vol._nrows, vol._ncols, vol._nslices) {

   double *source, *target;
   _allocate(2*this->get_nelements());
   for (source=vol._volume, target=_volume;
        target < _volume+2*get_nelements();
        source++, target += 2){
      *target     = *source;
      *(target+1) = 0.0;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume destructor
//---------------------------------------------------------------------------

MRI_Complex_Volume::~MRI_Complex_Volume() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Complex_Volume::fill_with(double fill_value){

   double *ptr;
   for(ptr=_volume; ptr<_volume+(2*get_nelements()); ptr++){
      *ptr = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::fill_real_with
// Fills the real-part of the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Complex_Volume::fill_real_with(double fill_value){
   
   double *ptr;
   for(ptr=_volume; ptr<_volume+(2*get_nelements()); ptr+=2){
      *ptr = fill_value;
   }

}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::fill_imag_with
// Fill the imaginary-part of the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Complex_Volume::fill_imag_with(double fill_value){

   double *ptr;
   for(ptr=_volume+1; ptr<_volume+(2*get_nelements()); ptr+=2){
      *ptr = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::real
// Returns the real part of the matrix.
//---------------------------------------------------------------------------

MRI_Double_Volume MRI_Complex_Volume::real(void) const {
   MRI_Double_Volume vol(_nrows,_ncols,_nslices);
   double *source, *target;
   for(source=_volume, target=vol._volume;
       source<_volume+(2*get_nelements());
       source+=2, target++){
      *target = *source;
   }
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::imag
// Returns the imaginary part of the matrix.
//---------------------------------------------------------------------------

MRI_Double_Volume MRI_Complex_Volume::imag(void) const {
   MRI_Double_Volume vol(_nrows,_ncols,_nslices);
   double *source, *target;
   for(source=_volume+1, target=vol._volume;
       source<_volume+(2*get_nelements());
       source+=2, target++){
      *target = *source;
   }
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Double_Volume::abs  
// Returns the modulus of the image.
//---------------------------------------------------------------------------

MRI_Double_Volume MRI_Complex_Volume::abs(void) const {
   MRI_Double_Volume vol(_nrows, _ncols, _nslices);
   double *source, *target;
   for(source=_volume, target=vol._volume;
       source<_volume+(2*get_nelements());
       source+=2, target++){
      *target = hypot(*source, *(source+1));
   }
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::angle
// Returns the angle (-PI..PI) of the matrix. 
//---------------------------------------------------------------------------

MRI_Double_Volume MRI_Complex_Volume::angle(void) const {
   MRI_Double_Volume vol(_nrows, _ncols, _nslices);
   double *source, *target;
   for(source=_volume, target=vol._volume;
       source<_volume+(2*get_nelements());
       source+=2, target++){
      *target = atan2(*(source+1),*source);
   }
   return vol;
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::operator=
// Assignment operator.
//---------------------------------------------------------------------------

MRI_Complex_Volume& MRI_Complex_Volume::operator=(const MRI_Complex_Volume& vol){
   if (this != &vol){
      _nrows = vol._nrows;
      _ncols = vol._ncols;
      _nslices = vol._nslices;
      (void)memcpy(_volume, vol._volume, this->size_in_bytes());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::operator*=
// Scales a matrix by (double) scalar. 
//---------------------------------------------------------------------------

MRI_Complex_Volume& MRI_Complex_Volume::operator*=(double a){
   double *ptr;
   for(ptr=_volume; ptr<_volume+2*get_nelements(); ptr++){
      *ptr *= a;
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::FFT2()
// Returns the 2D-FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Volume::FFT2(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned long sizes[2];
   unsigned int islice;
   double *slice;

   sizes[0] = _nrows;
   sizes[1] = _ncols;
   for(islice=0; islice<_nslices; islice++){
      slice = &(_volume[islice*_nrows*_ncols]);
      fourn(slice-1, sizes-1, 2, -1);
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::iFFT2
// Returns the 2D-inverse FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Volume::iFFT2(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned long sizes[2];
   unsigned int islice;
   double *slice;

   sizes[0] = _nrows;
   sizes[1] = _ncols;
   for(islice=0; islice<_nslices; islice++){
      slice = &(_volume[islice*_nrows*_ncols]);
      fourn(slice-1, sizes-1, 2, 1);
   }
   *this *= (1.0/(double)get_slicesize());
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::FFT3()
// Returns the 3D-FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Volume::FFT3(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned long sizes[3];
   sizes[0] = _nrows;
   sizes[1] = _ncols;
   sizes[2] = _nslices;
   fourn(_volume-1, sizes-1, 3, -1);
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::iFFT3
// Returns the 3D-inverse FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Volume::iFFT3(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned long sizes[3];
   sizes[0] = _nrows;
   sizes[1] = _ncols;
   sizes[2] = _nslices;
   fourn(_volume-1, sizes-1, 3, 1);
   *this *= (1.0/(double)get_nelements());
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::fftshift2
// Swaps first and fourth, second and third quadrants to move the
// zeroth lag to the centre of the spectrum.
//---------------------------------------------------------------------------

void MRI_Complex_Volume::fftshift2(void){
   unsigned int islice;
   double *slice;
   for(islice=0; islice<get_nslices(); islice++){
      slice = &(_volume[islice*_nrows*_ncols]);
      _fftshift((void *)slice, _nrows, _ncols, this->element_size_in_bytes());
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::fftshift3
// Swaps first and fourth, second and third quadrants to move the
// zeroth lag to the centre of the spectrum.
//---------------------------------------------------------------------------

void MRI_Complex_Volume::fftshift3(void){
// Swap octants
//   _fftshift((void *)_volume, _nrows, _ncols, this->element_size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Complex_Volume::display 
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_Complex_Volume::display(ostream& stream){
   unsigned int row, col, slice;
   for(slice=0; slice<_nslices; slice++){ 
      stream << endl;
      for(row=0; row<_nrows; row++){
         for(col=0; col<_ncols; col++){
            stream << "(" << setprecision(4) << this->real(row,col,slice) 
                   << "," << setprecision(4) << this->imag(row,col,slice) 
                   << ") ";
         }
      }
      stream << endl;
   }
}


