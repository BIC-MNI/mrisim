//===========================================================================
// MRIMATRIX.CXX
//
// R.Kwan
// August 31, 1995
//
// (C) Copyright 1995 by R.Kwan
//===========================================================================

//===========================================================================
// $Header: /private-cvsroot/simulation/mrisim/src/minc/mrimatrix.cxx,v 1.2 2004-08-10 15:24:19 bert Exp $
// $Log: mrimatrix.cxx,v $
// Revision 1.2  2004-08-10 15:24:19  bert
// Use atan2 rather than fatan2
//
// Revision 1.1  2003/05/30 16:43:09  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 3.1  1996/07/19  15:48:37  rkwan
// Release 3.1 update.
//
// Revision 2.5  1996/05/29  16:24:58  rkwan
// Release 2.5
//
// Revision 4.1  1995/12/22  20:18:38  rkwan
// Update for float type.
//
// Revision 3.1  1995/12/11  14:17:41  rkwan
// Updated for fast_iso_model.
//
//===========================================================================

#include "mrimatrix.h"

//===========================================================================
// MRI_Matrix
//===========================================================================

//---------------------------------------------------------------------------
// MRI_Matrix constructors
//---------------------------------------------------------------------------

MRI_Matrix::MRI_Matrix() {
   _nrows = 0;
   _ncols = 0;
}

MRI_Matrix::MRI_Matrix(unsigned int nrows, unsigned int ncols) {
   _nrows = nrows;
   _ncols = ncols;
}

//---------------------------------------------------------------------------
// MRI_Matrix destructor
//---------------------------------------------------------------------------

MRI_Matrix::~MRI_Matrix() {}

//---------------------------------------------------------------------------
// MRI_Matrix::is_same_size_as
// Returns TRUE if the two matrices have the same dimensions.
//---------------------------------------------------------------------------

int MRI_Matrix::is_same_size_as(const MRI_Matrix& mat) const {
   return (((_nrows == mat._nrows) && (_ncols == mat._ncols)) ? TRUE : FALSE);
}

//---------------------------------------------------------------------------
// MRI_Matrix::is_power_of_two
// Returns TRUE if the matrix dimensions are a power of two.
//---------------------------------------------------------------------------

int MRI_Matrix::is_power_of_two(void) const {
   return ((_power_of_two(_nrows) && _power_of_two(_ncols)) ? TRUE : FALSE);
}

//============================================================================
// MRI_Byte_Matrix
//============================================================================

//---------------------------------------------------------------------------
// MRI_Byte_Matrix constructors
//---------------------------------------------------------------------------

MRI_Byte_Matrix::MRI_Byte_Matrix() : MRI_Matrix() {
   _allocate(this->get_nelements());
}

MRI_Byte_Matrix::MRI_Byte_Matrix(unsigned int nrows, 
                             unsigned int ncols,
                             unsigned char fill_value) 
   : MRI_Matrix(nrows, ncols) {

   _allocate(this->get_nelements());
   this->fill_with(fill_value);
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix copy constructor 
//---------------------------------------------------------------------------

MRI_Byte_Matrix::MRI_Byte_Matrix(const MRI_Byte_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   _allocate(this->get_nelements());
   (void)memcpy(_matrix, mat._matrix, this->size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix type cast constructor
//---------------------------------------------------------------------------

MRI_Byte_Matrix::MRI_Byte_Matrix(const MRI_Short_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {
    int truncate = 0;
   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(this->get_nelements());

   for (n=0; n<len; n++){
      if (mat._matrix[n] < 0){
         _matrix[n] = 0;
         truncate   = 1;
      } else if (mat._matrix[n] > 255) {
         _matrix[n] = 255;
         truncate   = 1;
      } else {
         _matrix[n] = mat._matrix[n];
      }
   }

#ifdef DEBUG
   if (truncate == 1) {
      cerr << "WARNING:  Elements truncated in MRI_Byte_Matrix constructor."
           << endl;
   }
#endif

}

MRI_Byte_Matrix::MRI_Byte_Matrix(const MRI_Float_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);
  
   for (n=0; n<len; n++){
      _matrix[n] = (unsigned char) rint(mat._matrix[n]);
   } 
}

MRI_Byte_Matrix::MRI_Byte_Matrix(const MRI_Double_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);
  
   for (n=0; n<len; n++){
      _matrix[n] = (unsigned char) rint(mat._matrix[n]);
   } 
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix destructor
//---------------------------------------------------------------------------

MRI_Byte_Matrix::~MRI_Byte_Matrix() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Byte_Matrix::fill_with(unsigned char fill_value){
   (void)memset(_matrix, fill_value, this->get_nelements());
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix:operator=
// Assignment operator
//---------------------------------------------------------------------------

MRI_Byte_Matrix& MRI_Byte_Matrix::operator=(const MRI_Byte_Matrix& mat){
   if (this != &mat){
      _nrows = mat._nrows;
      _ncols = mat._ncols;
      (void)memcpy(_matrix, mat._matrix, this->get_nelements());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::operator*=
// Scales a matrix by a (double) scalar.
//---------------------------------------------------------------------------

MRI_Byte_Matrix& MRI_Byte_Matrix::operator*=(double a) { 
   unsigned int n;
   unsigned int len = this->get_nelements();
  
   for(n=0; n<len; n++){
      _matrix[n] = (unsigned char)rint(_matrix[n] * a);
   } 
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::operator*=
// Multiplies a matrix element-by-element by another matrix
//---------------------------------------------------------------------------

MRI_Byte_Matrix& MRI_Byte_Matrix::operator*=(const MRI_Byte_Matrix& x) { 
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
  
   for(n=0; n<len; n++){
      _matrix[n] = _matrix[n] * x._matrix[n];
   } 
   return *this;
}
//---------------------------------------------------------------------------
// MRI_Byte_Matrix::operator+=
// Add a matrix to the current matrix.
//---------------------------------------------------------------------------

MRI_Byte_Matrix& MRI_Byte_Matrix::operator+=(const MRI_Byte_Matrix& x){
   unsigned int n;
   unsigned int len = this->get_nelements();
 
#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
 
   for(n=0; n<len; n++){
      _matrix[n] = _matrix[n] + x._matrix[n];
   } 
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::saxpy
// Compute the scalar a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_Byte_Matrix& MRI_Byte_Matrix::saxpy(double a, const MRI_Byte_Matrix& x,
                                    const MRI_Byte_Matrix& y){
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
   assert(this->is_same_size_as(y));
#endif

   for(n=0; n<len; n++){
      _matrix[n] = (unsigned char)rint(a * x._matrix[n] + y._matrix[n]);
   }
   return *this;   
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::gaxpy
// Compute the general a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_Byte_Matrix& MRI_Byte_Matrix::gaxpy(const MRI_Byte_Matrix& a, 
                                    const MRI_Byte_Matrix& x,
                                    const MRI_Byte_Matrix& y){
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(a));
   assert(this->is_same_size_as(x));
   assert(this->is_same_size_as(y));
#endif

   for(n=0; n<len; n++){
      _matrix[n] = a._matrix[n] * x._matrix[n] + y._matrix[n];
   }
   return *this;   
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::operator !=
// Element-wise comparison.
//---------------------------------------------------------------------------

int MRI_Byte_Matrix::operator != (const MRI_Byte_Matrix& mat) const {
   unsigned int n = 0;
   unsigned int len = this->get_nelements();
   int equal = 1;
  
   do {
      if (_matrix[n] != mat._matrix[n]) equal = 0;
      n++;
   } while ((equal == 1) && (n<len));

   return equal;
}  

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::maximum 
// Returns the maximum element in the matrix.
//---------------------------------------------------------------------------

unsigned char MRI_Byte_Matrix::maximum(unsigned int row1, 
                                     unsigned int row2,
                                     unsigned int col1,
                                     unsigned int col2) const{

#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   unsigned char max = 0;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         if ((*this)(row,col) > max) max = (*this)(row,col);
      }
   }
   return max;
}
   
unsigned char MRI_Byte_Matrix::maximum(void) const {
   unsigned char max;
   unsigned int n;
   unsigned int len = this->get_nelements();
   
   for(max=0, n=0; n<len; n++){
      if (_matrix[n] > max) max = _matrix[n];
   }

   return max;
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::minimum
// Returns the minimum element in the matrix.
//---------------------------------------------------------------------------

unsigned char MRI_Byte_Matrix::minimum(unsigned int row1,
                                     unsigned int row2,
                                     unsigned int col1,
                                     unsigned int col2) const{

#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   unsigned char min = 255;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         if ((*this)(row,col) < min) min = (*this)(row,col);
      }
   }
   return min;
}
   
unsigned char MRI_Byte_Matrix::minimum(void) const{
   unsigned char min;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(min=255, n=0; n<len; n++){
      if (_matrix[n] < min) min = _matrix[n];
   }
   return min;
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::sum
// Returns the sum of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Byte_Matrix::sum(unsigned int row1,
                          unsigned int row2,
                          unsigned int col1,
                          unsigned int col2) const {

#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double sum = 0.0;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         sum += (double)((*this)(row,col));
      }
   }
   return sum;
}

double MRI_Byte_Matrix::sum(void) const {
   double sum = 0.0;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for (n=0; n<len; n++){
      sum += (double)_matrix[n];
   }
   return sum;
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::norm
// Returns the norm of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Byte_Matrix::norm(unsigned int row1,
                           unsigned int row2,
                           unsigned int col1,
                           unsigned int col2) const {

#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double avg = this->mean(row1,row2,col1,col2);
   double sum = 0.0;
   unsigned int row, col;

   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         sum += SQR((double)((*this)(row,col))-avg);
      }
   }
   return sqrt(sum);
}

double MRI_Byte_Matrix::norm(void) const {
   double avg = this->mean();
   double sum = 0.0;
   unsigned int n;
   unsigned len = this->get_nelements();

   for(n=0; n<len; n++){
      sum += SQR((double)_matrix[n]-avg);
   }
   return sqrt(sum);
}
 
//---------------------------------------------------------------------------
// MRI_Byte_Matrix::mean
// Returns the mean of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Byte_Matrix::mean(unsigned int row1,
                           unsigned int row2,
                           unsigned int col1,
                           unsigned int col2) const {

#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double sum = this->sum(row1,row2,col1,col2);
   return sum/(double)((row2-row1+1)*(col2-col1+1));
}

double MRI_Byte_Matrix::mean(void) const{
   double sum = this->sum(); 
   return sum/(double)get_nelements();
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::std
// Returns the standard deviation of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Byte_Matrix::std(unsigned int row1,
                           unsigned int row2,
                           unsigned int col1,
                           unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double norm = this->norm(row1,row2,col1,col2);
   return norm/sqrt((row2-row1+1)*(col2-col1+1)-1);
}

double MRI_Byte_Matrix::std(void) const{
   double norm = this->norm();
   return norm/sqrt(get_nelements()-1);
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::set_submatrix
// Copies a matrix into a submatrix
//---------------------------------------------------------------------------

void MRI_Byte_Matrix::set_submatrix(MRI_Byte_Matrix& mat, unsigned int row,
                                  unsigned int col) {
 
   unsigned char *source, *target;
   int offset = row*_ncols + col;

   if (this != &mat){
      for(source=mat._matrix, target=_matrix+offset;
          target<_matrix+get_nelements();
          source += mat._ncols, target += _ncols){
         memcpy(target, source, mat._ncols);
      }
   }
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::operator () 
// Returns a submatrix of the matrix.
//---------------------------------------------------------------------------

MRI_Byte_Matrix MRI_Byte_Matrix::operator() (unsigned int row1, unsigned int row2,
                                         unsigned int col1, unsigned int col2) const {

#ifdef DEBUG
   assert(row2 >= row1); 
   assert(col2 >= col1);
#endif

   unsigned int row, col;
   MRI_Byte_Matrix mat(row2-row1+1,col2-col1+1);
   for(row=row1; row<=row2; row++){
      for(col=col1; col<=col2; col++){
         mat(row-row1,col-col1) = (*this)(row,col);
      }
   }
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Byte_Matrix::display
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_Byte_Matrix::display(ostream& stream) const {
   unsigned int row, col;
   for(row=0; row<_nrows; row++){
      for(col=0; col<_ncols; col++){
         stream << (int)(*this)(row,col) << " ";
      }
      stream << endl;
   }
}

//============================================================================
// MRI_Short_Matrix
//============================================================================

//---------------------------------------------------------------------------
// MRI_Short_Matrix constructors
//---------------------------------------------------------------------------

MRI_Short_Matrix::MRI_Short_Matrix() : MRI_Matrix() {
   _allocate(this->get_nelements());
}

MRI_Short_Matrix::MRI_Short_Matrix(unsigned int nrows, 
                             unsigned int ncols,
                             short fill_value) 
   : MRI_Matrix(nrows, ncols) {

   _allocate(this->get_nelements());
   this->fill_with(fill_value);
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix copy constructor 
//---------------------------------------------------------------------------

MRI_Short_Matrix::MRI_Short_Matrix(const MRI_Short_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   _allocate(this->get_nelements());
   (void)memcpy(_matrix, mat._matrix, this->size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix type cast constructor
//---------------------------------------------------------------------------

MRI_Short_Matrix::MRI_Short_Matrix(const MRI_Byte_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);
   
   for(n=0; n<len; n++){
      _matrix[n] = (short) rint(mat._matrix[n]);
   }
}

MRI_Short_Matrix::MRI_Short_Matrix(const MRI_Float_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);

   for(n=0; n<len; n++){
      _matrix[n] = (short) rint(mat._matrix[n]);
   }
}

MRI_Short_Matrix::MRI_Short_Matrix(const MRI_Double_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);

   for(n=0; n<len; n++){
      _matrix[n] = (short) rint(mat._matrix[n]);
   }
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix destructor
//---------------------------------------------------------------------------

MRI_Short_Matrix::~MRI_Short_Matrix() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Short_Matrix::fill_with(short fill_value){
   unsigned int n;
   unsigned int len = this->get_nelements();
   for(n=0; n<len; n++){
      _matrix[n] = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix:operator=
// Assignment operator
//---------------------------------------------------------------------------

MRI_Short_Matrix& MRI_Short_Matrix::operator=(const MRI_Short_Matrix& mat){
   if (this != &mat){
      _nrows = mat._nrows;
      _ncols = mat._ncols;
      (void)memcpy(_matrix, mat._matrix, this->get_nelements());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::operator*=
// Scales a matrix by a (double) scalar.
//---------------------------------------------------------------------------

MRI_Short_Matrix& MRI_Short_Matrix::operator*=(double a) { 
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<len; n++){
      _matrix[n] = (short)rint(_matrix[n]*a);
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::operator*=
// Multiplies a matrix element-by-element by another matrix
//---------------------------------------------------------------------------

MRI_Short_Matrix& MRI_Short_Matrix::operator*=(const MRI_Short_Matrix& x) { 
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
  
   for(n=0; n<len; n++){
      _matrix[n] = _matrix[n] * x._matrix[n];
   } 
   return *this;
}
//---------------------------------------------------------------------------
// MRI_Short_Matrix::operator+=
// Add a matrix to the current matrix.
//---------------------------------------------------------------------------

MRI_Short_Matrix& MRI_Short_Matrix::operator+=(const MRI_Short_Matrix& x){
   unsigned int n;
   unsigned int len = this->get_nelements();
 
#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
 
   for(n=0; n<len; n++){
      _matrix[n] = _matrix[n] + x._matrix[n];
   } 
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::saxpy
// Computes the scalar a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_Short_Matrix& MRI_Short_Matrix::saxpy(double a, const MRI_Short_Matrix& x,
                                      const MRI_Short_Matrix& y){
   unsigned n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   // Check that all matrix operands have the same number of elements
   assert(len == x.get_nelements());
   assert(len == y.get_nelements());  
#endif
 
   for(n=0; n<len; n++){
      _matrix[n] = (short)rint(a * x._matrix[n] + y._matrix[n]);
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::gaxpy
// Computes the general a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_Short_Matrix& MRI_Short_Matrix::gaxpy(const MRI_Short_Matrix& a, 
                                      const MRI_Short_Matrix& x,
                                      const MRI_Short_Matrix& y){
   unsigned n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(a));
   assert(this->is_same_size_as(x));
   assert(this->is_same_size_as(y));
#endif
 
   for(n=0; n<len; n++){
      _matrix[n] = a._matrix[n] * x._matrix[n] + y._matrix[n];
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::operator !=
// Element-wise comparison.
//---------------------------------------------------------------------------

int MRI_Short_Matrix::operator != (const MRI_Short_Matrix& mat) const {
   unsigned int n = 0;
   unsigned int len = this->get_nelements();
   int equal = 1;

   do {
      if (_matrix[n] != mat._matrix[n]) equal = 0;
      n++;
   } while ((equal == 1) && (n<len));
   return equal;
}  

//---------------------------------------------------------------------------
// MRI_Short_Matrix::maximum 
// Returns the maximum element in the matrix.
//---------------------------------------------------------------------------

short MRI_Short_Matrix::maximum(unsigned int row1,
                                       unsigned int row2,
                                       unsigned int col1,
                                       unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   short max = SHRT_MIN;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         if ((*this)(row,col) > max) max = (*this)(row,col);
      }
   }
   return max;
}

short MRI_Short_Matrix::maximum(void) const {
   short max;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(max=SHRT_MIN, n=0; n<len; n++){
      if (_matrix[n] > max) max = _matrix[n];
   }
   return max;
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::minimum
// Returns the minimum element in the matrix.
//---------------------------------------------------------------------------

short MRI_Short_Matrix::minimum(unsigned int row1,
                                       unsigned int row2,
                                       unsigned int col1,
                                       unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   short min = SHRT_MAX;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         if ((*this)(row,col) < min) min = (*this)(row,col);
      }
   }
   return min;
}
   
short MRI_Short_Matrix::minimum(void) const {
   short min;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(min=SHRT_MAX, n=0; n<len; n++){
      if (_matrix[n] < min) min = _matrix[n];
   }
   return min;
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::sum
// Returns the sum of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Short_Matrix::sum(unsigned int row1,
                           unsigned int row2,
                           unsigned int col1,
                           unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double sum = 0.0;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         sum += (double)((*this)(row,col));
      }
   }
   return sum;
}

double MRI_Short_Matrix::sum(void) const {
   double sum = 0.0;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<len; n++){
      sum += (double)_matrix[n];
   }
   return sum;
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::norm
// Returns the norm of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Short_Matrix::norm(unsigned int row1,
                            unsigned int row2,
                            unsigned int col1,
                            unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double avg = this->mean(row1,row2,col1,col2);
   double sum = 0.0;
   unsigned int row, col;

   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         sum += SQR((double)((*this)(row,col))-avg);
      }
   }
   return sqrt(sum);
}

double MRI_Short_Matrix::norm(void) const {
   double avg = this->mean();
   double sum = 0.0;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<len; n++){
      sum += SQR((double)_matrix[n]-avg);
   }
   return sqrt(sum);
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::mean
// Returns the mean of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Short_Matrix::mean(unsigned int row1,
                            unsigned int row2,
                            unsigned int col1,
                            unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double sum = this->sum(row1,row2,col1,col2);
   return sum/(double)((row2-row1+1)*(col2-col1+1));
}

double MRI_Short_Matrix::mean(void) const{
   double sum = this->sum(); 
   return sum/(double)get_nelements();
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::std
// Returns the standard deviation of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Short_Matrix::std(unsigned int row1,
                            unsigned int row2,
                            unsigned int col1,
                            unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double norm = this->norm(row1,row2,col1,col2);
   return norm/sqrt((row2-row1+1)*(col2-col1+1)-1);
}

double MRI_Short_Matrix::std(void) const{
   double norm = this->norm();
   return norm/sqrt(get_nelements()-1);
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::set_submatrix
// Copies a matrix into a submatrix
//---------------------------------------------------------------------------

void MRI_Short_Matrix::set_submatrix(MRI_Short_Matrix& mat, unsigned int row,
                                  unsigned int col) {
  
   short *source, *target;
   int offset = row*_ncols + col;

   if (this != &mat){
      for(source=mat._matrix, target=_matrix+offset;
          target<_matrix+get_nelements();
          source += mat._ncols, target += _ncols){
         memcpy(target, source, mat._ncols);
      }
   }
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::operator () 
// Returns a submatrix of the matrix.
//---------------------------------------------------------------------------

MRI_Short_Matrix MRI_Short_Matrix::operator() (unsigned int row1, unsigned int row2,
                                         unsigned int col1, unsigned int col2) const {

#ifdef DEBUG
   assert(row2 >= row1); 
   assert(col2 >= col1);
#endif

   unsigned int row, col;
   MRI_Short_Matrix mat(row2-row1+1,col2-col1+1);
   for(row=row1; row<=row2; row++){
      for(col=col1; col<=col2; col++){
         mat(row-row1,col-col1) = (*this)(row,col);
      }
   }
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Short_Matrix::display
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_Short_Matrix::display(ostream& stream) const {
   unsigned int row, col;
   for(row=0; row<_nrows; row++){
      for(col=0; col<_ncols; col++){
         stream << (short int)(*this)(row,col) << " ";
      }
      stream << endl;
   }
}
 
//===========================================================================
// MRI_Float_Matrix
//===========================================================================

//---------------------------------------------------------------------------
// MRI_Float_Matrix constructors
//---------------------------------------------------------------------------

MRI_Float_Matrix::MRI_Float_Matrix() : MRI_Matrix() {
   _allocate(this->get_nelements());
}

MRI_Float_Matrix::MRI_Float_Matrix(unsigned int nrows, unsigned int ncols,
                                 float fill) :
    MRI_Matrix(nrows, ncols) {

   _allocate(this->get_nelements());
   this->fill_with(fill);
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix copy constructor
//---------------------------------------------------------------------------

MRI_Float_Matrix::MRI_Float_Matrix(const MRI_Float_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   _allocate(this->get_nelements());
   (void)memcpy(_matrix, mat._matrix, this->size_in_bytes());

}

//---------------------------------------------------------------------------
// MRI_Float_Matrix type cast constructor
//---------------------------------------------------------------------------

MRI_Float_Matrix::MRI_Float_Matrix(const MRI_Byte_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);

   for(n=0; n<len; n++){
      _matrix[n] = (float)mat._matrix[n];
   }
}

MRI_Float_Matrix::MRI_Float_Matrix(const MRI_Short_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);

   for(n=0; n<len; n++){
      _matrix[n] = (float)mat._matrix[n];
   }
} 

MRI_Float_Matrix::MRI_Float_Matrix(const MRI_Double_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);

   for(n=0; n<len; n++){
      _matrix[n] = (float)mat._matrix[n];
   }
} 

//---------------------------------------------------------------------------
// MRI_Float_Matrix destructor
//---------------------------------------------------------------------------

MRI_Float_Matrix::~MRI_Float_Matrix() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Float_Matrix::fill_with(float fill_value) {
   unsigned int n;
   unsigned int len = this->get_nelements();
   for(n=0; n<len; n++){
      _matrix[n] = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::operator=
// Assignment operator.
//---------------------------------------------------------------------------

MRI_Float_Matrix& MRI_Float_Matrix::operator=(const MRI_Float_Matrix& mat){
   if (this != &mat){
      _nrows = mat._nrows;
      _ncols = mat._ncols;
      (void)memcpy(_matrix, mat._matrix, this->size_in_bytes());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::operator*=
// Scales a matrix by a (float) scalar.
//---------------------------------------------------------------------------

MRI_Float_Matrix& MRI_Float_Matrix::operator*=(double a){
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<len; n++){
      _matrix[n] *= a;
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::operator*=
// Multiplies a matrix element-by-element by another matrix
//---------------------------------------------------------------------------

MRI_Float_Matrix& MRI_Float_Matrix::operator*=(const MRI_Float_Matrix& x) { 
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
  
   for(n=0; n<len; n++){
      _matrix[n] *= x._matrix[n];
   } 
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::operator/=
// Divides a matrix element-by-element by another matrix
//---------------------------------------------------------------------------

MRI_Float_Matrix& MRI_Float_Matrix::operator/=(const MRI_Float_Matrix& x) {
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif

   for(n=0; n<len; n++){
      _matrix[n] /= x._matrix[n];
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::operator+=
// Add a matrix to the current matrix.
//---------------------------------------------------------------------------

MRI_Float_Matrix& MRI_Float_Matrix::operator+=(const MRI_Float_Matrix& x){
   unsigned int n;
   unsigned int len = this->get_nelements();
 
#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
 
   for(n=0; n<len; n++){
      _matrix[n] += x._matrix[n];
   } 
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::saxpy
// Compute the scalar a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_Float_Matrix& MRI_Float_Matrix::saxpy(double a, const MRI_Float_Matrix& x,
                                      const MRI_Float_Matrix& y){
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   // Check that all matrix operands have the same number of elements
   assert(len == x.get_nelements());
   assert(len == y.get_nelements());
#endif

   for(n=0; n<len; n++){
      _matrix[n] = a * x._matrix[n] + y._matrix[n];
   }
   return *this;   
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::gaxpy
// Compute the general a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_Float_Matrix& MRI_Float_Matrix::gaxpy(const MRI_Float_Matrix& a, 
                                      const MRI_Float_Matrix& x,
                                      const MRI_Float_Matrix& y){
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(a));
   assert(this->is_same_size_as(x));
   assert(this->is_same_size_as(y));
#endif

   for(n=0; n<len; n++){
      _matrix[n] = a._matrix[n] * x._matrix[n] + y._matrix[n];
   }
   return *this;   
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::operator!=
// Element-wise comparison.
//---------------------------------------------------------------------------

int MRI_Float_Matrix::operator != (const MRI_Float_Matrix& mat) const {
   unsigned int n = 0;
   unsigned int len = this->get_nelements();
   int equal = 1;

   do {
      if (_matrix[n] != mat._matrix[n]) equal = 0;
      n++;
   } while ((equal == 1) && (n<len));
   return equal;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::maximum
// Returns the maximum element in the matrix.
//---------------------------------------------------------------------------

double MRI_Float_Matrix::maximum(unsigned int row1,
                                unsigned int row2,
                                unsigned int col1,
                                unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double max = FLT_MIN;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         if ((*this)(row,col) > max) max = (*this)(row,col);
      }
   }
   return max;
}

double MRI_Float_Matrix::maximum(void) const {
   double max;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(max=FLT_MIN, n=0; n<len; n++){
      if (_matrix[n] > max) max = _matrix[n];
   }
   return max;
}
   
//---------------------------------------------------------------------------
// MRI_Float_Matrix::minimum
// Returns the minimum element in the matrix.
//---------------------------------------------------------------------------

double MRI_Float_Matrix::minimum(unsigned int row1,
                                unsigned int row2,
                                unsigned int col1,
                                unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double min = FLT_MAX;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         if ((*this)(row,col) < min) min = (*this)(row,col);
      }
   }
   return min;
}

double MRI_Float_Matrix::minimum(void) const{
   double min;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(min=FLT_MAX, n=0; n<len; n++){
      if (_matrix[n] < min) min = _matrix[n];
   }
   return min;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::sum
// Returns the sum of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Float_Matrix::sum(unsigned int row1,
                            unsigned int row2,
                            unsigned int col1,
                            unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double sum = 0.0;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         sum += (double)((*this)(row,col));
      }
   }
   return sum;
}

double MRI_Float_Matrix::sum(void) const {
   unsigned int n;
   unsigned int len = this->get_nelements();
   double sum = 0.0;
   for (n=0; n<len; n++){
      sum += _matrix[n];
   }
   return sum;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::norm
// Returns the norm of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Float_Matrix::norm(unsigned int row1,
                             unsigned int row2,
                             unsigned int col1,
                             unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double avg = this->mean(row1,row2,col1,col2);
   double sum = 0.0;
   unsigned int row, col;

   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         sum += SQR((double)((*this)(row,col))-avg);
      }
   }
   return sqrt(sum);
}

double MRI_Float_Matrix::norm(void) const {
   double avg = this->mean();
   double sum = 0.0;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<len; n++){
      sum += SQR(_matrix[n]-avg);
   }
   return sqrt(sum);
}
 
//---------------------------------------------------------------------------
// MRI_Float_Matrix::mean
// Returns the mean of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Float_Matrix::mean(unsigned int row1,
                             unsigned int row2,
                             unsigned int col1,
                             unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double sum = this->sum(row1,row2,col1,col2);
   return sum/((row2-row1+1)*(col2-col1+1));
}

double MRI_Float_Matrix::mean(void) const{
   double sum = this->sum(); 
   return sum/get_nelements();
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::std
// Returns the standard deviation of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Float_Matrix::std(unsigned int row1,
                            unsigned int row2,
                            unsigned int col1,
                            unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double norm = this->norm(row1,row2,col1,col2);
   return norm/sqrt((row2-row1+1)*(col2-col1+1)-1);
}

double MRI_Float_Matrix::std(void) const{
   double norm = this->norm();
   return norm/sqrt(get_nelements()-1);
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::FFT
// Returns the 1D-FFT of the rows of the matrix.
//---------------------------------------------------------------------------

MRI_FComplex_Matrix MRI_Float_Matrix::FFT(void) const {
   MRI_FComplex_Matrix mat(*this);
   mat.FFT();
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::iFFT
// Returns the 1D-inverse FFT of the rows of the matrix.
//---------------------------------------------------------------------------

MRI_FComplex_Matrix MRI_Float_Matrix::iFFT(void) const {
   MRI_FComplex_Matrix mat(*this);
   mat.iFFT();
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::FFT2
// Returns the 2D-FFT of the matrix.
//---------------------------------------------------------------------------

MRI_FComplex_Matrix MRI_Float_Matrix::FFT2(void) const {
   MRI_FComplex_Matrix mat(*this); 
   mat.FFT2();
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::iFFT2
// Returns the 2D-inverse FFT of the matrix.
//---------------------------------------------------------------------------

MRI_FComplex_Matrix MRI_Float_Matrix::iFFT2(void) const {
   MRI_FComplex_Matrix mat(*this);
   mat.iFFT2();
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::fftshift
// Swaps first and fourth, second and third quadrants to move the
// zeroth lag to the centre of the spectrum.
//---------------------------------------------------------------------------

void MRI_Float_Matrix::fftshift(void){
   _fftshift((void *)_matrix, _nrows, _ncols, this->element_size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::set_submatrix
// Copies a matrix into a submatrix.
//---------------------------------------------------------------------------

void MRI_Float_Matrix::set_submatrix(MRI_Float_Matrix& mat, unsigned int row,
                                    unsigned int col){

   float *source, *target;
   int offset = row*_ncols+col;

   if (this != &mat){
      for(source=mat._matrix, target=_matrix+offset;
          target<_matrix+get_nelements();
          source += mat._ncols, target += _ncols){
         (void)memcpy(target, source, mat._ncols*sizeof(float));
      }
   }
} 

//---------------------------------------------------------------------------
// MRI_Float_Matrix::operator()
// Returns a submatrix of the matrix.
//---------------------------------------------------------------------------

MRI_Float_Matrix MRI_Float_Matrix::operator() (unsigned int row1, 
                                           unsigned int row2,
                                           unsigned int col1, 
                                           unsigned int col2) const {

#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   unsigned int row, col;
   MRI_Float_Matrix mat(row2-row1+1,col2-col1+1);
   for(row=row1; row<=row2; row++){
      for(col=col1; col<=col2; col++){
         mat(row-row1,col-col1) = (*this)(row,col);
      }
   }
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Float_Matrix::display
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_Float_Matrix::display(ostream& stream) const {
   unsigned int row, col;
   for(row=0; row<_nrows; row++){
      for(col=0; col<_ncols; col++){
         stream << setprecision(4) << (*this)(row,col) << " ";
      }
      stream << endl;
   }
}

//===========================================================================
// MRI_Double_Matrix
//===========================================================================

//---------------------------------------------------------------------------
// MRI_Double_Matrix constructors
//---------------------------------------------------------------------------

MRI_Double_Matrix::MRI_Double_Matrix() : MRI_Matrix() {
   _allocate(this->get_nelements());
}

MRI_Double_Matrix::MRI_Double_Matrix(unsigned int nrows, unsigned int ncols,
                                 double fill) :
    MRI_Matrix(nrows, ncols) {

   _allocate(this->get_nelements());
   this->fill_with(fill);
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix copy constructor
//---------------------------------------------------------------------------

MRI_Double_Matrix::MRI_Double_Matrix(const MRI_Double_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   _allocate(this->get_nelements());
   (void)memcpy(_matrix, mat._matrix, this->size_in_bytes());

}

//---------------------------------------------------------------------------
// MRI_Double_Matrix type cast constructor
//---------------------------------------------------------------------------

MRI_Double_Matrix::MRI_Double_Matrix(const MRI_Byte_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);

   for(n=0; n<len; n++){
      _matrix[n] = (double)mat._matrix[n];
   }
}

MRI_Double_Matrix::MRI_Double_Matrix(const MRI_Short_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);

   for(n=0; n<len; n++){
      _matrix[n] = (double)mat._matrix[n];
   }
} 

MRI_Double_Matrix::MRI_Double_Matrix(const MRI_Float_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();

   _allocate(len);

   for(n=0; n<len; n++){
      _matrix[n] = (double)mat._matrix[n];
   }
} 

//---------------------------------------------------------------------------
// MRI_Double_Matrix destructor
//---------------------------------------------------------------------------

MRI_Double_Matrix::~MRI_Double_Matrix() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Double_Matrix::fill_with(double fill_value) {
   unsigned int n;
   unsigned int len = this->get_nelements();
   for(n=0; n<len; n++){
      _matrix[n] = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::operator=
// Assignment operator.
//---------------------------------------------------------------------------

MRI_Double_Matrix& MRI_Double_Matrix::operator=(const MRI_Double_Matrix& mat){
   if (this != &mat){
      _nrows = mat._nrows;
      _ncols = mat._ncols;
      (void)memcpy(_matrix, mat._matrix, this->size_in_bytes());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::operator*=
// Scales a matrix by a (double) scalar.
//---------------------------------------------------------------------------

MRI_Double_Matrix& MRI_Double_Matrix::operator*=(double a){
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<len; n++){
      _matrix[n] *= a;
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::operator*=
// Multiplies a matrix element-by-element by another matrix
//---------------------------------------------------------------------------

MRI_Double_Matrix& MRI_Double_Matrix::operator*=(const MRI_Double_Matrix& x) { 
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
  
   for(n=0; n<len; n++){
      _matrix[n] *= x._matrix[n];
   } 
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::operator/=
// Divides a matrix element-by-element by another matrix
//---------------------------------------------------------------------------

MRI_Double_Matrix& MRI_Double_Matrix::operator/=(const MRI_Double_Matrix& x) {
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
 
   for(n=0; n<len; n++){
      _matrix[n] /= x._matrix[n];
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::operator+=
// Add a matrix to the current matrix.
//---------------------------------------------------------------------------

MRI_Double_Matrix& MRI_Double_Matrix::operator+=(const MRI_Double_Matrix& x){
   unsigned int n;
   unsigned int len = this->get_nelements();
 
#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
 
   for(n=0; n<len; n++){
      _matrix[n] = _matrix[n] + x._matrix[n];
   } 
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::saxpy
// Compute the scalar a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_Double_Matrix& MRI_Double_Matrix::saxpy(double a, const MRI_Double_Matrix& x,
                                        const MRI_Double_Matrix& y){
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   // Check that all matrix operands have the same number of elements
   assert(len == x.get_nelements());
   assert(len == y.get_nelements());
#endif

   for(n=0; n<len; n++){
      _matrix[n] = a * x._matrix[n] + y._matrix[n];
   }
   return *this;   
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::gaxpy
// Compute the general a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_Double_Matrix& MRI_Double_Matrix::gaxpy(const MRI_Double_Matrix& a, 
                                        const MRI_Double_Matrix& x,
                                        const MRI_Double_Matrix& y){
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(a));
   assert(this->is_same_size_as(x));
   assert(this->is_same_size_as(y));
#endif

   for(n=0; n<len; n++){
      _matrix[n] = a._matrix[n] * x._matrix[n] + y._matrix[n];
   }
   return *this;   
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::operator!=
// Element-wise comparison.
//---------------------------------------------------------------------------

int MRI_Double_Matrix::operator != (const MRI_Double_Matrix& mat) const {
   unsigned int n = 0;
   unsigned int len = this->get_nelements();
   int equal = 1;

   do {
      if (_matrix[n] != mat._matrix[n]) equal = 0;
      n++;
   } while ((equal == 1) && (n<len));
   return equal;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::maximum
// Returns the maximum element in the matrix.
//---------------------------------------------------------------------------

double MRI_Double_Matrix::maximum(unsigned int row1,
                                unsigned int row2,
                                unsigned int col1,
                                unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double max = DBL_MIN;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         if ((*this)(row,col) > max) max = (*this)(row,col);
      }
   }
   return max;
}

double MRI_Double_Matrix::maximum(void) const {
   double max;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(max=DBL_MIN, n=0; n<len; n++){
      if (_matrix[n] > max) max = _matrix[n];
   }
   return max;
}
   
//---------------------------------------------------------------------------
// MRI_Double_Matrix::minimum
// Returns the minimum element in the matrix.
//---------------------------------------------------------------------------

double MRI_Double_Matrix::minimum(unsigned int row1,
                                unsigned int row2,
                                unsigned int col1,
                                unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double min = DBL_MAX;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         if ((*this)(row,col) < min) min = (*this)(row,col);
      }
   }
   return min;
}

double MRI_Double_Matrix::minimum(void) const{
   double min;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(min=DBL_MAX, n=0; n<len; n++){
      if (_matrix[n] < min) min = _matrix[n];
   }
   return min;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::sum
// Returns the sum of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Double_Matrix::sum(unsigned int row1,
                            unsigned int row2,
                            unsigned int col1,
                            unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double sum = 0.0;
   unsigned int row, col;
   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         sum += (double)((*this)(row,col));
      }
   }
   return sum;
}

double MRI_Double_Matrix::sum(void) const {
   unsigned int n;
   unsigned int len = this->get_nelements();
   double sum = 0.0;
   for (n=0; n<len; n++){
      sum += _matrix[n];
   }
   return sum;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::norm
// Returns the norm of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Double_Matrix::norm(unsigned int row1,
                             unsigned int row2,
                             unsigned int col1,
                             unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double avg = this->mean(row1,row2,col1,col2);
   double sum = 0.0;
   unsigned int row, col;

   for (row=row1; row<=row2; row++){
      for (col=col1; col<=col2; col++){
         sum += SQR((double)((*this)(row,col))-avg);
      }
   }
   return sqrt(sum);
}

double MRI_Double_Matrix::norm(void) const {
   double avg = this->mean();
   double sum = 0.0;
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<len; n++){
      sum += SQR(_matrix[n]-avg);
   }
   return sqrt(sum);
}
 
//---------------------------------------------------------------------------
// MRI_Double_Matrix::mean
// Returns the mean of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Double_Matrix::mean(unsigned int row1,
                             unsigned int row2,
                             unsigned int col1,
                             unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double sum = this->sum(row1,row2,col1,col2);
   return sum/((row2-row1+1)*(col2-col1+1));
}

double MRI_Double_Matrix::mean(void) const{
   double sum = this->sum(); 
   return sum/get_nelements();
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::std
// Returns the standard deviation of the matrix elements.
//---------------------------------------------------------------------------

double MRI_Double_Matrix::std(unsigned int row1,
                            unsigned int row2,
                            unsigned int col1,
                            unsigned int col2) const {
#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   double norm = this->norm(row1,row2,col1,col2);
   return norm/sqrt((row2-row1+1)*(col2-col1+1)-1);
}

double MRI_Double_Matrix::std(void) const{
   double norm = this->norm();
   return norm/sqrt(get_nelements()-1);
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::FFT
// Returns the 1D-FFT of the rows of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Matrix MRI_Double_Matrix::FFT(void) const {
   MRI_Complex_Matrix mat(*this);
   mat.FFT();
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::iFFT
// Returns the 1D-inverse FFT of the rows of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Matrix MRI_Double_Matrix::iFFT(void) const {
   MRI_Complex_Matrix mat(*this);
   mat.iFFT();
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::FFT2
// Returns the 2D-FFT of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Matrix MRI_Double_Matrix::FFT2(void) const {
   MRI_Complex_Matrix mat(*this); 
   mat.FFT2();
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::iFFT2
// Returns the 2D-inverse FFT of the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Matrix MRI_Double_Matrix::iFFT2(void) const {
   MRI_Complex_Matrix mat(*this);
   mat.iFFT2();
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::fftshift
// Swaps first and fourth, second and third quadrants to move the
// zeroth lag to the centre of the spectrum.
//---------------------------------------------------------------------------

void MRI_Double_Matrix::fftshift(void){
   _fftshift((void *)_matrix, _nrows, _ncols, this->element_size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::set_submatrix
// Copies a matrix into a submatrix.
//---------------------------------------------------------------------------

void MRI_Double_Matrix::set_submatrix(MRI_Double_Matrix& mat, unsigned int row,
                                    unsigned int col){

   double *source, *target;
   int offset = row*_ncols+col;

   if (this != &mat){
      for(source=mat._matrix, target=_matrix+offset;
          target<_matrix+get_nelements();
          source += mat._ncols, target += _ncols){
         (void)memcpy(target, source, mat._ncols*sizeof(double));
      }
   }
} 

//---------------------------------------------------------------------------
// MRI_Double_Matrix::operator()
// Returns a submatrix of the matrix.
//---------------------------------------------------------------------------

MRI_Double_Matrix MRI_Double_Matrix::operator() (unsigned int row1, unsigned int row2,
                                             unsigned int col1, unsigned int col2) const {

#ifdef DEBUG
   assert(row2 >= row1);
   assert(col2 >= col1);
#endif

   unsigned int row, col;
   MRI_Double_Matrix mat(row2-row1+1,col2-col1+1);
   for(row=row1; row<=row2; row++){
      for(col=col1; col<=col2; col++){
         mat(row-row1,col-col1) = (*this)(row,col);
      }
   }
   return mat;
}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::display
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_Double_Matrix::display(ostream& stream) const {
   unsigned int row, col;
   for(row=0; row<_nrows; row++){
      for(col=0; col<_ncols; col++){
         stream << setprecision(4) << (*this)(row,col) << " ";
      }
      stream << endl;
   }
}

//===========================================================================
// MRI_FComplex_Matrix
//===========================================================================

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix constructors
//---------------------------------------------------------------------------

MRI_FComplex_Matrix::MRI_FComplex_Matrix() : MRI_Matrix() {
   _allocate(2*this->get_nelements());
}

MRI_FComplex_Matrix::MRI_FComplex_Matrix(unsigned int nrows, 
                                         unsigned int ncols,
                                         float fill_value) :
    MRI_Matrix(nrows, ncols) {

   _allocate(2*this->get_nelements());
   this->fill_with(fill_value);
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix copy constructor
//---------------------------------------------------------------------------

MRI_FComplex_Matrix::MRI_FComplex_Matrix(const MRI_FComplex_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   _allocate(2*this->get_nelements());
   (void)memcpy(_matrix, mat._matrix, this->size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix type cast constructor
//---------------------------------------------------------------------------

MRI_FComplex_Matrix::MRI_FComplex_Matrix(const MRI_Byte_Matrix& mat) :
   MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();
   
   _allocate(2*len);

   for(n=0; n<len; n++){
      _matrix[2*n]   = (float)mat._matrix[n];
      _matrix[2*n+1] = 0.0;
   }
}

MRI_FComplex_Matrix::MRI_FComplex_Matrix(const MRI_Float_Matrix& mat) :
   MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();
   
   _allocate(2*len);

   for(n=0; n<len; n++){
      _matrix[2*n]   = mat._matrix[n];
      _matrix[2*n+1] = 0.0;
   }
}

MRI_FComplex_Matrix::MRI_FComplex_Matrix(const MRI_Double_Matrix& mat) :
   MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();
   
   _allocate(2*len);

   for(n=0; n<len; n++){
      _matrix[2*n]   = (float)mat._matrix[n];
      _matrix[2*n+1] = 0.0;
   }
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix destructor
//---------------------------------------------------------------------------

MRI_FComplex_Matrix::~MRI_FComplex_Matrix() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::fill_with(float fill_value){
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<2*len; n++){
      _matrix[n] = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::fill_real_with
// Fills the real-part of the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::fill_real_with(float fill_value){
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<2*len; n+=2){
      _matrix[n] = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::fill_imag_with
// Fill the imaginary-part of the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::fill_imag_with(float fill_value){
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=1; n<2*len; n+=2){
      _matrix[n] = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::real
// Returns the real part of the matrix.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::real(MRI_Float_Matrix& mat) const {
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->get_nrows() == mat.get_nrows());
   assert(this->get_ncols() == mat.get_ncols());
#endif

   for(n=0; n<len; n++){
      mat._matrix[n] = _matrix[2*n];
   }

}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::imag
// Returns the imaginary part of the matrix.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::imag(MRI_Float_Matrix& mat) const {
   unsigned int n;
   unsigned int len = this->get_nelements();
 
#ifdef DEBUG
   assert(this->get_nrows() == mat.get_nrows());
   assert(this->get_ncols() == mat.get_ncols());
#endif

   for(n=0; n<len; n++){
      mat._matrix[n] = _matrix[2*n+1];
   }

}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::abs  
// Returns the modulus of the image.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::abs(MRI_Float_Matrix& mat) const {
   unsigned int n;
   unsigned int len = this->get_nelements();
   
#ifdef DEBUG
   assert(this->get_nrows() == mat.get_nrows());
   assert(this->get_ncols() == mat.get_ncols());
#endif

   for(n=0; n<len; n++){
      mat._matrix[n] = hypotf(_matrix[2*n], _matrix[2*n+1]);
   }

}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::angle
// Returns the angle (-PI..PI) of the matrix. 
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::angle(MRI_Float_Matrix& mat) const {
   unsigned int n;
   unsigned int len = this->get_nelements();
  
#ifdef DEBUG
   assert(this->get_nrows() == mat.get_nrows());
   assert(this->get_ncols() == mat.get_ncols());
#endif

   for(n=0; n<len; n++){
      mat._matrix[n] = (float) atan2(_matrix[2*n+1], _matrix[2*n]);
   }

}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::get_real_min_max
// Gets the minimum and maximum real values in the matrix.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::get_real_min_max(double &min, double &max) const {

   double tmp1, tmp2, swap;
   min = FLT_MAX; max = FLT_MIN;

   unsigned int n, len = this->get_nelements();

   for (n=0; n<2*len; n+=4){
      tmp1 = _matrix[n];
      tmp2 = _matrix[n+2];
      if (tmp2 < tmp1){
         swap = tmp1; tmp1 = tmp2; tmp2 = swap;
      }
      if (tmp1 < min) min = tmp1;
      if (tmp2 > max) max = tmp2;
   }
   if (!(len % 2)) { // if number of elements is odd
      tmp1 = _matrix[2*len-2];
      if (tmp1 < min) min = tmp1;
      if (tmp1 > max) max = tmp1;
   }
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::get_imag_min_max
// Gets the minimum and maximum imaginary values in the matrix.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::get_imag_min_max(double &min, double &max) const {

   double tmp1, tmp2, swap;
   min = FLT_MAX; max = FLT_MIN;

   unsigned int n, len = this->get_nelements();

   for (n=1; n<2*len; n+=4){
      tmp1 = _matrix[n];
      tmp2 = _matrix[n+2];
      if (tmp2 < tmp1){
         swap = tmp1; tmp1 = tmp2; tmp2 = swap;
      }
      if (tmp1 < min) min = tmp1;
      if (tmp2 > max) max = tmp2;
   }
   if (!(len % 2)) { // if number of elements is odd
      tmp1 = _matrix[2*len-2];
      if (tmp1 < min) min = tmp1;
      if (tmp1 > max) max = tmp1;
   }
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::get_abs_min_max
// Gets the minimum and maximum magnitude of the matrix.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::get_abs_min_max(double &min, double &max) const {

   double tmp1, tmp2, swap;
   min = FLT_MAX; max = FLT_MIN;

   unsigned int n, len = this->get_nelements();

   for (n=0; n<2*len; n+=4){
      tmp1 = hypotf(_matrix[n],   _matrix[n+1]);
      tmp2 = hypotf(_matrix[n+2], _matrix[n+3]);
      if (tmp2 < tmp1){
         swap = tmp1; tmp1 = tmp2; tmp2 = swap;
      }
      if (tmp1 < min) min = tmp1;
      if (tmp2 > max) max = tmp2;
   }
   if (!(len % 2)) { // if number of elements is odd
      tmp1 = hypotf(_matrix[2*len-2], _matrix[2*len-1]);
      if (tmp1 < min) min = tmp1;
      if (tmp1 > max) max = tmp1;
   }
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::get_angle_min_max
// Gets the minimum and maximum phase of the matrix.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::get_angle_min_max(double &min, double &max) const {

   double tmp1, tmp2, swap;
   min = FLT_MAX; max = FLT_MIN;

   unsigned int n, len = this->get_nelements();

   for (n=0; n<2*len; n+=4){
      tmp1 = atan2(_matrix[n+1], _matrix[n]);
      tmp2 = atan2(_matrix[n+3], _matrix[n+2]);
      if (tmp2 < tmp1){
         swap = tmp1; tmp1 = tmp2; tmp2 = swap;
      }
      if (tmp1 < min) min = tmp1;
      if (tmp2 > max) max = tmp2;
   }
   if (!(len % 2)) { // if number of elements is odd
      tmp1 = atan2(_matrix[2*len-1], _matrix[2*len-2]);
      if (tmp1 < min) min = tmp1;
      if (tmp1 > max) max = tmp1;
   }
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::operator=
// Assignmnt operator.
//---------------------------------------------------------------------------

MRI_FComplex_Matrix& MRI_FComplex_Matrix::operator=(const MRI_FComplex_Matrix& mat){
   if (this != &mat){
      _nrows = mat._nrows;
      _ncols = mat._ncols;
      (void)memcpy(_matrix, mat._matrix, this->size_in_bytes());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::operator*=
// Scales a matrix by (double) scaler. 
//---------------------------------------------------------------------------

MRI_FComplex_Matrix& MRI_FComplex_Matrix::operator*=(double a){
   unsigned int n;
   unsigned int len = 2*this->get_nelements();
  
   for(n=0; n<len; n++){
      _matrix[n] *= a;
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::operator*=
// Multiplies a matrix element-by-element by another matrix
//---------------------------------------------------------------------------

MRI_FComplex_Matrix& MRI_FComplex_Matrix::operator*=(const MRI_FComplex_Matrix& x) { 
   unsigned int n;
   unsigned int len = 2*this->get_nelements();
   double   re, im;

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
  
   for(n=0; n<len; n+=2){
      re = _matrix[n]*x._matrix[n]   - _matrix[n+1]*x._matrix[n+1];
      im = _matrix[n]*x._matrix[n+1] + _matrix[n+1]*x._matrix[n];
      _matrix[n]   = (float)re;
      _matrix[n+1] = (float)im;
   } 
   return *this;
}

MRI_FComplex_Matrix& MRI_FComplex_Matrix::operator*=(const MRI_Float_Matrix& x) {
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif

   for(n=0; n<len; n++){
      _matrix[2*n]   *= x._matrix[n];
      _matrix[2*n+1] *= x._matrix[n];
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::operator/=
// Divides a matrix element-by-element by another matrix.
//---------------------------------------------------------------------------

MRI_FComplex_Matrix& MRI_FComplex_Matrix::operator/=(const MRI_Float_Matrix& x)
{
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif

   for(n=0; n<len; n++){
      _matrix[2*n]   /= x._matrix[n];
      _matrix[2*n+1] /= x._matrix[n];
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::operator+=
// Add a matrix to the current matrix.
//---------------------------------------------------------------------------

MRI_FComplex_Matrix& MRI_FComplex_Matrix::operator+=(const MRI_FComplex_Matrix& x){
   unsigned int n;
   unsigned int len = 2*this->get_nelements();
 
#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
 
   for(n=0; n<len; n++){
      _matrix[n] = _matrix[n] + x._matrix[n];
   } 
   return *this;
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::saxpy
// Compute the scalar a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_FComplex_Matrix& MRI_FComplex_Matrix::saxpy(double areal, double aimag, 
                                          const MRI_FComplex_Matrix& x,
                                          const MRI_FComplex_Matrix& y){
   unsigned int n;
   unsigned int len = 2*this->get_nelements();
   double   re, im;

#ifdef DEBUG
   assert(this->is_same_size_as(x));
   assert(this->is_same_size_as(y));
#endif

   // Use temporary re, im in case x or y point to this.
   for(n=0; n<len; n+=2){
      re = areal*x._matrix[n]   - aimag*x._matrix[n+1] + y._matrix[n];
      im = areal*x._matrix[n+1] + aimag*x._matrix[n]   + y._matrix[n+1];
      _matrix[n] = (float)re;
      _matrix[n] = (float)im;
   }
   return *this;   
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::gaxpy
// Compute the general a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_FComplex_Matrix& MRI_FComplex_Matrix::gaxpy(const MRI_FComplex_Matrix& a,
                                          const MRI_FComplex_Matrix& x,
                                          const MRI_FComplex_Matrix& y){
   unsigned int n;
   unsigned int len = 2*this->get_nelements();
   double   re, im;

#ifdef DEBUG
   assert(this->is_same_size_as(x));
   assert(this->is_same_size_as(y));
#endif

   // Use temporary re, im in case a, x or y point to this.
   for(n=0; n<len; n+=2){
      re = a._matrix[n]*x._matrix[n]   - a._matrix[n+1]*x._matrix[n+1] + 
           y._matrix[n];
      im = a._matrix[n]*x._matrix[n+1] + a._matrix[n+1]*x._matrix[n]   + 
           y._matrix[n+1];
      _matrix[n]   = (float)re;
      _matrix[n+1] = (float)im;
   }
   return *this;   
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::FFT
// Compute the row-wise 1D FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::FFT(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned int  irow;
   float        *row;

   for(irow=0; irow<_nrows; irow++){
      row = &(_matrix[irow*_ncols]);
      four1f(row-1, _ncols, -1);
   }
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::iFFT
// Compute the row-wise 1D FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::iFFT(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned int  irow;
   float        *row;

   for(irow=0; irow<_nrows; irow++){
      row = &(_matrix[irow*_ncols]);
      four1f(row-1, _ncols, 1);
   }
   *this *= (1.0/(float)get_ncols());
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::FFT2()
// Returns the 2D-FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::FFT2(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned long sizes[2];
   sizes[0] = _nrows;
   sizes[1] = _ncols;
   fournf(_matrix-1, sizes-1, 2, -1);
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::iFFT2
// Returns the 2D-inverse FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::iFFT2(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned long sizes[2];
   sizes[0] = _nrows;
   sizes[1] = _ncols;
   fournf(_matrix-1, sizes-1, 2, 1);
   *this *= (1.0/(float)get_nelements());
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::fftshift
// Swaps first and fourth, second and third quadrants to move the
// zeroth lag to the centre of the spectrum.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::fftshift(void){
   _fftshift((void *)_matrix, _nrows, _ncols, this->element_size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::set_submatrix
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::set_submatrix(MRI_FComplex_Matrix& mat, 
                                     unsigned int row,
                                     unsigned int col) {

   float *source, *target;
   int offset = 2*(row*_ncols + col);
   
   if (this != &mat){
      for(source=mat._matrix, target=_matrix+offset;
          target<_matrix+get_nelements();
           source += mat._ncols, target += _ncols){
         memcpy(target, source, 2*mat._ncols*sizeof(float));
      }
   }
}
                             
//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::set_submatrix
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::set_submatrix(MRI_Float_Matrix& mat, 
                                     unsigned int row,
                                     unsigned int col) {

   unsigned int i,j,m,n;
   for(i=row, m=0; m<mat.get_nrows(); i++, m++){
      for(j=col, n=0; n<mat.get_ncols(); j++, n++){
         this->real(i,j) = (float)mat(m,n);
         this->imag(i,j) = 0.0;
      }
   }
}

//---------------------------------------------------------------------------
// MRI_FComplex_Matrix::display 
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_FComplex_Matrix::display(ostream& stream) const {
   unsigned int row, col;
   for(row=0; row<_nrows; row++){
      for(col=0; col<_ncols; col++){
         stream << "(" << setprecision(4) << this->real(row,col) 
                << "," << setprecision(4) << this->imag(row,col) 
                << ") ";
      }
      stream << endl;
   }
}

//===========================================================================
// MRI_Complex_Matrix
//===========================================================================

//---------------------------------------------------------------------------
// MRI_Complex_Matrix constructors
//---------------------------------------------------------------------------

MRI_Complex_Matrix::MRI_Complex_Matrix() : MRI_Matrix() {
   _allocate(2*this->get_nelements());
}

MRI_Complex_Matrix::MRI_Complex_Matrix(unsigned int nrows, unsigned int ncols,
                                   double fill_value) :
    MRI_Matrix(nrows, ncols) {

   _allocate(2*this->get_nelements());
   this->fill_with(fill_value);
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix copy constructor
//---------------------------------------------------------------------------

MRI_Complex_Matrix::MRI_Complex_Matrix(const MRI_Complex_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   _allocate(2*mat.get_nelements());
   (void)memcpy(_matrix, mat._matrix, this->size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix type cast constructor
//---------------------------------------------------------------------------

MRI_Complex_Matrix::MRI_Complex_Matrix(const MRI_FComplex_Matrix& mat) :
    MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = mat.get_nelements();

   _allocate(2*len);

   for(n=0; n<2*len; n++){
      _matrix[n] = (double)mat._matrix[n];
   }
}

MRI_Complex_Matrix::MRI_Complex_Matrix(const MRI_Byte_Matrix& mat) :
   MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();
   
   _allocate(2*len);

   for(n=0; n<len; n++){
      _matrix[2*n]   = (double)mat._matrix[n];
      _matrix[2*n+1] = 0.0;
   }
}

MRI_Complex_Matrix::MRI_Complex_Matrix(const MRI_Float_Matrix& mat) :
   MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();
   
   _allocate(2*len);

   for(n=0; n<len; n++){
      _matrix[2*n]   = (double)mat._matrix[n];
      _matrix[2*n+1] = 0.0;
   }
}

MRI_Complex_Matrix::MRI_Complex_Matrix(const MRI_Double_Matrix& mat) :
   MRI_Matrix(mat._nrows, mat._ncols) {

   unsigned int n;
   unsigned int len = this->get_nelements();
   
   _allocate(2*len);

   for(n=0; n<len; n++){
      _matrix[2*n]   = mat._matrix[n];
      _matrix[2*n+1] = 0.0;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix destructor
//---------------------------------------------------------------------------

MRI_Complex_Matrix::~MRI_Complex_Matrix() {
   _deallocate();
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::fill_with
// Fills the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::fill_with(double fill_value){
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<2*len; n++){
      _matrix[n] = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::fill_real_with
// Fills the real-part of the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::fill_real_with(double fill_value){
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=0; n<2*len; n+=2){
      _matrix[n] = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::fill_imag_with
// Fill the imaginary-part of the matrix with the given fill value.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::fill_imag_with(double fill_value){
   unsigned int n;
   unsigned int len = this->get_nelements();

   for(n=1; n<2*len; n+=2){
      _matrix[n] = fill_value;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::real
// Returns the real part of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::real(MRI_Double_Matrix& mat) const {
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->get_nrows() == mat.get_nrows());
   assert(this->get_ncols() == mat.get_ncols());
#endif

   for(n=0; n<len; n++){
      mat._matrix[n] = _matrix[2*n];
   }

}

//---------------------------------------------------------------------------
// MRI_Double_Matrix::imag
// Returns the imaginary part of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::imag(MRI_Double_Matrix& mat) const {
   unsigned int n;
   unsigned int len = this->get_nelements();
 
#ifdef DEBUG
   assert(this->get_nrows() == mat.get_nrows());
   assert(this->get_ncols() == mat.get_ncols());
#endif

   for(n=0; n<len; n++){
      mat._matrix[n] = _matrix[2*n+1];
   }

}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::abs  
// Returns the modulus of the image.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::abs(MRI_Double_Matrix& mat) const {
   unsigned int n;
   unsigned int len = this->get_nelements();
  
#ifdef DEBUG
   assert(this->get_nrows() == mat.get_nrows());
   assert(this->get_ncols() == mat.get_ncols());
#endif
 
   for(n=0; n<len; n++){
      mat._matrix[n] = hypot(_matrix[2*n], _matrix[2*n+1]);
   }

}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::angle
// Returns the angle (-PI..PI) of the matrix. 
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::angle(MRI_Double_Matrix& mat) const {
   unsigned int n;
   unsigned int len = this->get_nelements();
  
#ifdef DEBUG
   assert(this->get_nrows() == mat.get_nrows());
   assert(this->get_ncols() == mat.get_ncols());
#endif

   for(n=0; n<len; n++){
      mat._matrix[n] = atan2(_matrix[2*n+1], _matrix[2*n]);
   }

}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::get_real_min_max
// Gets the minimum and maximum real values in the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::get_real_min_max(double &min, double &max) const {

   double tmp1, tmp2, swap;
   min = DBL_MAX; max = DBL_MIN;

   unsigned int n, len = this->get_nelements();

   for (n=0; n<2*len; n+=4){
      tmp1 = _matrix[n];
      tmp2 = _matrix[n+2];
      if (tmp2 < tmp1){
         swap = tmp1; tmp1 = tmp2; tmp2 = swap;
      }
      if (tmp1 < min) min = tmp1;
      if (tmp2 > max) max = tmp2;
   }
   if (!(len % 2)) { // if number of elements is odd
      tmp1 = _matrix[2*len-2];
      if (tmp1 < min) min = tmp1;
      if (tmp1 > max) max = tmp1;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::get_imag_min_max
// Gets the minimum and maximum imaginary values in the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::get_imag_min_max(double &min, double &max) const {

   double tmp1, tmp2, swap;
   min = DBL_MAX; max = DBL_MIN;

   unsigned int n, len = this->get_nelements();

   for (n=1; n<2*len; n+=4){
      tmp1 = _matrix[n];
      tmp2 = _matrix[n+2];
      if (tmp2 < tmp1){
         swap = tmp1; tmp1 = tmp2; tmp2 = swap;
      }
      if (tmp1 < min) min = tmp1;
      if (tmp2 > max) max = tmp2;
   }
   if (!(len % 2)) { // if number of elements is odd
      tmp1 = _matrix[2*len-2];
      if (tmp1 < min) min = tmp1;
      if (tmp1 > max) max = tmp1;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::get_abs_min_max
// Gets the minimum and maximum magnitude of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::get_abs_min_max(double &min, double &max) const {

   double tmp1, tmp2, swap;
   min = DBL_MAX; max = DBL_MIN;

   unsigned int n, len = this->get_nelements();

   for (n=0; n<2*len; n+=4){
      tmp1 = hypot(_matrix[n],   _matrix[n+1]);
      tmp2 = hypot(_matrix[n+2], _matrix[n+3]);
      if (tmp2 < tmp1){
         swap = tmp1; tmp1 = tmp2; tmp2 = swap;
      }
      if (tmp1 < min) min = tmp1;
      if (tmp2 > max) max = tmp2;
   }
   if (!(len % 2)) { // if number of elements is odd
      tmp1 = hypot(_matrix[2*len-2], _matrix[2*len-1]);
      if (tmp1 < min) min = tmp1;
      if (tmp1 > max) max = tmp1;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::get_angle_min_max
// Gets the minimum and maximum phase of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::get_angle_min_max(double &min, double &max) const {

   double tmp1, tmp2, swap;
   min = DBL_MAX; max = DBL_MIN;

   unsigned int n, len = this->get_nelements();

   for (n=0; n<2*len; n+=4){
      tmp1 = atan2(_matrix[n+1], _matrix[n]);
      tmp2 = atan2(_matrix[n+3], _matrix[n+2]);
      if (tmp2 < tmp1){
         swap = tmp1; tmp1 = tmp2; tmp2 = swap;
      }
      if (tmp1 < min) min = tmp1;
      if (tmp2 > max) max = tmp2;
   }
   if (!(len % 2)) { // if number of elements is odd
      tmp1 = atan2(_matrix[2*len-1], _matrix[2*len-2]);
      if (tmp1 < min) min = tmp1;
      if (tmp1 > max) max = tmp1;
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::operator=
// Assignmnt operator.
//---------------------------------------------------------------------------

MRI_Complex_Matrix& MRI_Complex_Matrix::operator=(const MRI_Complex_Matrix& mat){
   if (this != &mat){
      _nrows = mat._nrows;
      _ncols = mat._ncols;
      (void)memcpy(_matrix, mat._matrix, this->size_in_bytes());
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::operator*=
// Scales a matrix by (double) scaler. 
//---------------------------------------------------------------------------

MRI_Complex_Matrix& MRI_Complex_Matrix::operator*=(double a){
   unsigned int n;
   unsigned int len = 2*this->get_nelements();
  
   for(n=0; n<len; n++){
      _matrix[n] *= a;
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::operator*=
// Multiplies a matrix element-by-element by another matrix
//---------------------------------------------------------------------------

MRI_Complex_Matrix& MRI_Complex_Matrix::operator*=(const MRI_Complex_Matrix& x) { 
   unsigned int n;
   unsigned int len = 2*this->get_nelements();
   double   re, im;

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
  
   for(n=0; n<len; n+=2){
      re = _matrix[n]*x._matrix[n]   - _matrix[n+1]*x._matrix[n+1];
      im = _matrix[n]*x._matrix[n+1] + _matrix[n+1]*x._matrix[n];
      _matrix[n]   = re;
      _matrix[n+1] = im;
   } 
   return *this;
}

MRI_Complex_Matrix& MRI_Complex_Matrix::operator*=(const MRI_Double_Matrix& x) {

   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif

   for(n=0; n<len; n++){
      _matrix[2*n]   *= x._matrix[n];
      _matrix[2*n+1] *= x._matrix[n];
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::operator/=
// Divides a matrix element-by-element by another matrix.
//---------------------------------------------------------------------------

MRI_Complex_Matrix& MRI_Complex_Matrix::operator/=(const MRI_Double_Matrix& x) {
   unsigned int n;
   unsigned int len = this->get_nelements();

#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif

   for(n=0; n<len; n++){
      _matrix[2*n]   /= x._matrix[n];
      _matrix[2*n+1] /= x._matrix[n];
   }
   return *this;
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::operator+=
// Add a matrix to the current matrix.
//---------------------------------------------------------------------------

MRI_Complex_Matrix& MRI_Complex_Matrix::operator+=(const MRI_Complex_Matrix& x){
   unsigned int n;
   unsigned int len = 2*this->get_nelements();
 
#ifdef DEBUG
   assert(this->is_same_size_as(x));
#endif
 
   for(n=0; n<len; n++){
      _matrix[n] = _matrix[n] + x._matrix[n];
   } 
   return *this;
}


//---------------------------------------------------------------------------
// MRI_Complex_Matrix::saxpy
// Compute the scalar a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Matrix& MRI_Complex_Matrix::saxpy(double areal, double aimag, 
                                          const MRI_Complex_Matrix& x,
                                          const MRI_Complex_Matrix& y){
   unsigned int n;
   unsigned int len = 2*this->get_nelements();
   double   re, im;

#ifdef DEBUG
   assert(this->is_same_size_as(x));
   assert(this->is_same_size_as(y));
#endif

   // Use temporary re, im in case x or y point to this.
   for(n=0; n<len; n+=2){
      re = areal*x._matrix[n]   - aimag*x._matrix[n+1] + y._matrix[n];
      im = areal*x._matrix[n+1] + aimag*x._matrix[n]   + y._matrix[n+1];
      _matrix[n] = re;
      _matrix[n] = im;
   }
   return *this;   
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::gaxpy
// Compute the general a * x + y results and saves it in the matrix.
//---------------------------------------------------------------------------

MRI_Complex_Matrix& MRI_Complex_Matrix::gaxpy(const MRI_Complex_Matrix& a,
                                          const MRI_Complex_Matrix& x,
                                          const MRI_Complex_Matrix& y){
   unsigned int n;
   unsigned int len = 2*this->get_nelements();
   double   re, im;

#ifdef DEBUG
   assert(this->is_same_size_as(x));
   assert(this->is_same_size_as(y));
#endif

   // Use temporary re, im in case a, x or y point to this.
   for(n=0; n<len; n+=2){
      re = a._matrix[n]*x._matrix[n]   - a._matrix[n+1]*x._matrix[n+1] + 
           y._matrix[n];
      im = a._matrix[n]*x._matrix[n+1] + a._matrix[n+1]*x._matrix[n]   + 
           y._matrix[n+1];
      _matrix[n]   = re;
      _matrix[n+1] = im;
   }
   return *this;   
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::FFT
// Compute the row-wise 1D FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::FFT(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned int  irow;
   double        *row;

   for(irow=0; irow<_nrows; irow++){
      row = &(_matrix[irow*_ncols]);
      four1(row-1, _ncols, -1);
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::iFFT
// Compute the row-wise 1D FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::iFFT(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned int  irow;
   double        *row;

   for(irow=0; irow<_nrows; irow++){
      row = &(_matrix[irow*_ncols]);
      four1(row-1, _ncols, 1);
   }
   *this *= (1.0/(double)get_ncols());
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::FFT2()
// Returns the 2D-FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::FFT2(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned long sizes[2];
   sizes[0] = _nrows;
   sizes[1] = _ncols;
   fourn(_matrix-1, sizes-1, 2, -1);
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::iFFT2
// Returns the 2D-inverse FFT of the matrix.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::iFFT2(void){

#ifdef DEBUG
   assert(this->is_power_of_two());
#endif

   unsigned long sizes[2];
   sizes[0] = _nrows;
   sizes[1] = _ncols;
   fourn(_matrix-1, sizes-1, 2, 1);
   *this *= (1.0/(double)get_nelements());
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::fftshift
// Swaps first and fourth, second and third quadrants to move the
// zeroth lag to the centre of the spectrum.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::fftshift(void){
   _fftshift((void *)_matrix, _nrows, _ncols, this->element_size_in_bytes());
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::set_submatrix
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::set_submatrix(MRI_Complex_Matrix& mat, unsigned int row,
                                     unsigned int col) {

   double *source, *target;
   int offset = 2*(row*_ncols + col);
   
   if (this != &mat){
      for(source=mat._matrix, target=_matrix+offset;
          target<_matrix+get_nelements();
           source += mat._ncols, target += _ncols){
         memcpy(target, source, 2*mat._ncols*sizeof(double));
      }
   }
}
                             
//---------------------------------------------------------------------------
// MRI_Complex_Matrix::set_submatrix
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::set_submatrix(MRI_Byte_Matrix& mat, unsigned int row,
                                     unsigned int col) {

   unsigned int i,j,m,n;
   for(i=row, m=0; m<mat.get_nrows(); i++, m++){
      for(j=col, n=0; n<mat.get_ncols(); j++, n++){
         this->real(i,j) = (double)mat(m,n);
         this->imag(i,j) = 0.0;
      }
   }
}

//---------------------------------------------------------------------------
// MRI_Complex_Matrix::display 
// Outputs the matrix elements to an output stream.
//---------------------------------------------------------------------------

void MRI_Complex_Matrix::display(ostream& stream) const {
   unsigned int row, col;
   for(row=0; row<_nrows; row++){
      for(col=0; col<_ncols; col++){
         stream << "(" << setprecision(4) << this->real(row,col) 
                << "," << setprecision(4) << this->imag(row,col) 
                << ") ";
      }
      stream << endl;
   }
}

