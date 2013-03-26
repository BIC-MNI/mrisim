#ifndef __VECTOR_H
#define __VECTOR_H

/*****************************************************************************
 * 
 * VECTOR.H
 *
 * VIO_Vector classes
 *
 * R. Kwan
 * August 8, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/vector.h,v 1.1 2003-05-30 16:43:14 bert Exp $
 * $Log: vector.h,v $
 * Revision 1.1  2003-05-30 16:43:14  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.2  1995/12/11  14:31:30  rkwan
 * Added arbitrary rotation axis.
 * Updated for fast_iso_model.
 *
 * Revision 2.1  1995/08/08  14:56:40  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include <mrisim/mrisim.h>

/*****************************************************************************
 * Types 
 *****************************************************************************/

enum Axis {X_AXIS=0, Y_AXIS=1, Z_AXIS=2};   // axis specifier for 3D vectors
class VIO_Vector;
extern double abs(VIO_Vector&);

/*****************************************************************************
 * VIO_Vector Class
 * Defines a generalized N-dimensional vector object and basic operations
 * on vectors.
 *****************************************************************************/

class VIO_Vector {
   public:
      VIO_Vector(int len);                       // constructors
      VIO_Vector(int len, double data[]);
      VIO_Vector(int len, double a);
      VIO_Vector(const VIO_Vector& v);               // copy constructor

      virtual ~VIO_Vector();                     // destructor

      double& operator[](int index);         // access indiv. vector elements
     
      void print(ostream& output) const;     // output functions

      friend double abs(VIO_Vector& v);
      friend int    dim(VIO_Vector& v);

   protected:
      double  *_vector;                      // storage for vector elements
      int     _length;                       // number of vector elements

};

/*****************************************************************************
 * Vector_3D Class
 * 3D Euclidean space vector object.  Inherits from VIO_Vector class.
 *****************************************************************************/

class Vector_3D : public VIO_Vector {
   public:
      Vector_3D();                           // constructors
      Vector_3D(double data[]);
      Vector_3D(double x, double y, double z);
      Vector_3D(double a);
      Vector_3D(const Vector_3D& v);

      virtual ~Vector_3D() {}                // destructor

      double& operator[](Axis axis);         // access vector components

      Vector_3D& operator=(const Vector_3D& v);  // assignment   
      Vector_3D& operator=(double a);
      Vector_3D& operator+=(const Vector_3D& v);
      Vector_3D& operator*=(const Vector_3D& v);
      Vector_3D& operator+=(double a);
      Vector_3D& operator*=(double a);

      void rotate_x(Degrees angle);       // 3D vector rotation
      void rotate_y(Degrees angle);
      void rotate_z(Degrees angle);
      void rotate(Axis axis, Degrees angle);
      void rotate(Vector_3D& axis, Degrees angle);   
};

#endif

