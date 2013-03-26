/*****************************************************************************
 * 
 * VECTOR.CXX
 *
 * VIO_Vector class member functions
 *
 * R. Kwan
 * August 8, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

// $Header: /private-cvsroot/simulation/mrisim/src/signal/vector.cxx,v 1.2 2008-11-06 10:58:23 rotor Exp $
// $Log: vector.cxx,v $
// Revision 1.2  2008-11-06 10:58:23  rotor
//  * fixed includes for iostream and friends
//  * updated for new release (1.0.2)
//
// Revision 1.1  2003/05/30 16:43:14  bert
// Initial checkin, mrisim 3.1 from Remi Kwan's home directory
//
// Revision 2.2  1995/12/11  14:31:49  rkwan
// Added arbitrary rotation axis.
// Updated for fast_iso_model.
//
// Revision 2.1  1995/08/08  14:56:47  rkwan
// Updated for new spin model and pulse sequence.
//

#include <math.h>
#include <iostream>
#include "vector.h"

/*****************************************************************************
 * VIO_Vector Class Methods
 *****************************************************************************/

/*****************************************************************************
 * VIO_Vector Constructor
 *****************************************************************************/


VIO_Vector::VIO_Vector(int len){
   int i;

   _length = len;
   _vector = new double[_length];

   for (i=0; i<_length; i++){
      _vector[i] = 0.0;
   }
}

VIO_Vector::VIO_Vector(int len, double data[]){
   int i;
   
   _length = len;
   _vector = new double[_length];

   for (i=0; i<_length; i++){
      _vector[i] = data[i];
   }
}

VIO_Vector::VIO_Vector(int len, double a){
   int i;
   
   _length = len;
   _vector = new double[_length];
   for (i=0; i<_length; i++){
      _vector[i] = a;
   }
}

/*****************************************************************************
 * VIO_Vector copy constructor
 *****************************************************************************/

VIO_Vector::VIO_Vector(const VIO_Vector& v){
   int i;

   _length = v._length;
   _vector = new double[_length];

   for (i=0; i<_length; i++){
      _vector[i] = v._vector[i];
   }
}

/*****************************************************************************
 * VIO_Vector Destructor
 *****************************************************************************/

VIO_Vector::~VIO_Vector(){
   delete[] _vector;
}

/*****************************************************************************
 * VIO_Vector::operator[]
 * Allows access to an element of the vector.
 *****************************************************************************/
  
double& VIO_Vector::operator[](int index){
   if (index > _length-1)
      return _vector[_length-1];
   else
      return _vector[index];
}

/*****************************************************************************
 * VIO_Vector::print
 * Prints a vector
 *****************************************************************************/

void VIO_Vector::print(ostream& output) const {
   int i;

   output << "[";
   for (i=0; i<_length-1; i++){
      output << _vector[i] << ", ";
   }
   output << _vector[_length-1] << "]" << endl;
}

/*****************************************************************************
 * VIO_Vector friend functions
 *****************************************************************************/

/*****************************************************************************
 * abs(VIO_Vector&)
 * Returns the absolute value of the vector
 *****************************************************************************/

double abs(VIO_Vector& v){
   int i;
   double modulus;

   for (modulus=0, i=0; i<v._length; i++){
      modulus += SQR(v[i]);
   }
   return sqrt(modulus);
}

/*****************************************************************************
 * dim(VIO_Vector&)
 * Returns the dimension of the vector
 *****************************************************************************/

int dim(VIO_Vector& v){
   return v._length;
}

/*****************************************************************************
 * Vector_3D Class Methods
 *****************************************************************************/

/*****************************************************************************
 * Vector_3D constructor
 *****************************************************************************/

Vector_3D::Vector_3D() : VIO_Vector(3) { }

Vector_3D::Vector_3D(double data[]) : VIO_Vector(3, data) { }

Vector_3D::Vector_3D(double x, double y, double z) : VIO_Vector(3){
   _vector[X_AXIS] = x;
   _vector[Y_AXIS] = y;
   _vector[Z_AXIS] = z;
}

Vector_3D::Vector_3D(double a) : VIO_Vector(3, a) { }

Vector_3D::Vector_3D(const Vector_3D& v) : VIO_Vector(3){
   _vector[X_AXIS] = v._vector[X_AXIS];
   _vector[Y_AXIS] = v._vector[Y_AXIS];
   _vector[Z_AXIS] = v._vector[Z_AXIS];
}

/*****************************************************************************
 * Vector_3D::operator[]
 * Allows access to an element of the vector.
 *****************************************************************************/

double& Vector_3D::operator[](Axis axis){
   return _vector[axis];
}

/*****************************************************************************
 * Vector_3D::operator=
 * Assignment operator.
 *****************************************************************************/

Vector_3D& Vector_3D::operator=(const Vector_3D& v){

   if (this != &v){
      _vector[X_AXIS] = v._vector[X_AXIS];
      _vector[Y_AXIS] = v._vector[Y_AXIS];
      _vector[Z_AXIS] = v._vector[Z_AXIS];
   }
   return *this;
}

Vector_3D& Vector_3D::operator=(double a){
   _vector[X_AXIS] = a;
   _vector[Y_AXIS] = a;
   _vector[Z_AXIS] = a;
   return *this;
}

/*****************************************************************************
 * Vector_3D::operator+=(const Vector_3D&)
 * Element-wise vector addition.
 *****************************************************************************/

Vector_3D& Vector_3D::operator+=(const Vector_3D& v){
   _vector[X_AXIS] += v._vector[X_AXIS];
   _vector[Y_AXIS] += v._vector[Y_AXIS];
   _vector[Z_AXIS] += v._vector[Z_AXIS];
   return *this;
}

/*****************************************************************************
 * Vector_3D::operator+=(double)
 * Add scalar to each vector element.
 *****************************************************************************/

Vector_3D& Vector_3D::operator+=(double a){
   _vector[X_AXIS] += a;
   _vector[Y_AXIS] += a;
   _vector[Z_AXIS] += a;
   return *this;
}

/*****************************************************************************
 * Vector_3D::operator*=(const Vector_3D&)
 * Element-wise vector multiplication
 *****************************************************************************/

Vector_3D& Vector_3D::operator*=(const Vector_3D& v){
   _vector[X_AXIS] *= v._vector[X_AXIS];
   _vector[Y_AXIS] *= v._vector[Y_AXIS];
   _vector[Z_AXIS] *= v._vector[Z_AXIS];
   return *this;
}

/*****************************************************************************
 * Vector_3D::operator*=(double)
 * Multiply scalar by each vector element.
 *****************************************************************************/

Vector_3D& Vector_3D::operator*=(double a){
   _vector[X_AXIS] *= a;
   _vector[Y_AXIS] *= a;
   _vector[Z_AXIS] *= a;
   return *this;
}

/*****************************************************************************
 * Vector_3D::rotate_x
 * Rotate a 3D vector about the x-axis.
 *****************************************************************************/

void Vector_3D::rotate_x(Degrees angle){
   double sina, cosa;
   double y, z;

   if (angle == 90){
      sina = 1;  cosa = 0;
   } else if (angle == 180){
      sina = 0;  cosa = -1;
   } else {
      sina = sin(DEG_TO_RAD(angle));
      cosa = cos(DEG_TO_RAD(angle));
   }
   y    = _vector[Y_AXIS];
   z    = _vector[Z_AXIS];

   _vector[Y_AXIS] =  y*cosa + z*sina;
   _vector[Z_AXIS] = -y*sina + z*cosa;
}

/*****************************************************************************
 * Vector_3D::rotate_y
 * Rotate a 3D vector about the y-axis.
 *****************************************************************************/

void Vector_3D::rotate_y(Degrees angle){
   double sina, cosa;
   double x, z;

   if (angle == 90){
      sina = 1;  cosa = 0;
   } else if (angle == 180){
      sina = 0;  cosa = -1;
   } else {
      sina = sin(DEG_TO_RAD(angle));
      cosa = cos(DEG_TO_RAD(angle));
   }
   x    = _vector[X_AXIS];
   z    = _vector[Z_AXIS];

   _vector[X_AXIS] =  x*cosa - z*sina;
   _vector[Z_AXIS] =  x*sina + z*cosa;
}

/*****************************************************************************
 * Vector_3D::rotate_z
 * Rotate a 3D vector about the z-axis.
 *****************************************************************************/

void Vector_3D::rotate_z(Degrees angle){
   double sina, cosa;
   double x, y;

   if (angle == 90){
      sina = 1;  cosa = 0;
   } else if (angle == 180){
      sina = 0;  cosa = -1;
   } else {
      sina = sin(DEG_TO_RAD(angle));
      cosa = cos(DEG_TO_RAD(angle));
   }

   x    = _vector[X_AXIS];
   y    = _vector[Y_AXIS];

   _vector[X_AXIS] =  x*cosa + y*sina;
   _vector[Y_AXIS] = -x*sina + y*cosa;
}

/*****************************************************************************
 * Vector_3D::rotate
 * Rotate a 3D vector about the specified axis.
 *****************************************************************************/

void Vector_3D::rotate(Axis axis, Degrees angle){

   switch(axis){
     case X_AXIS:
        rotate_x(angle);
        break;
     case Y_AXIS:
        rotate_y(angle);
        break;
     case Z_AXIS:
        rotate_z(angle);
        break;
  }
}

void Vector_3D::rotate(Vector_3D& axis, Degrees angle){

   double x = _vector[X_AXIS];
   double y = _vector[Y_AXIS];
   double z = _vector[Z_AXIS];

   double sina = sin(DEG_TO_RAD(angle));
   double uxsin = axis[X_AXIS]*sina;
   double uysin = axis[Y_AXIS]*sina;
   double uzsin = axis[Z_AXIS]*sina;

   double cosa = cos(DEG_TO_RAD(angle));
   double uxcos = (1-cosa)*axis[X_AXIS];
   double uycos = (1-cosa)*axis[Y_AXIS];
   double uzcos = (1-cosa)*axis[Z_AXIS];

   _vector[X_AXIS] = (axis[X_AXIS]*uxcos+cosa)*x
                   + (axis[Y_AXIS]*uxcos+uzsin)*y
                   + (axis[Z_AXIS]*uxcos-uysin)*z;
   _vector[Y_AXIS] = (axis[X_AXIS]*uycos-uzsin)*x
                   + (axis[Y_AXIS]*uycos+cosa)*y
                   + (axis[Z_AXIS]*uycos+uxsin)*z;
   _vector[Z_AXIS] = (axis[X_AXIS]*uzcos+uysin)*x
                   + (axis[Y_AXIS]*uzcos-uxsin)*y
                   + (axis[Z_AXIS]*uzcos+cosa)*z;
}
