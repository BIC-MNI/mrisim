#ifndef __MRISIM_H
#define __MRISIM_H

//==========================================================================
// MRISIM.H
// MRISIM global types and constants.
// Inherits from:
// Base class to: 
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/mrisim.h,v 1.2 2008-11-06 10:58:23 rotor Exp $
 * $Log: mrisim.h,v $
 * Revision 1.2  2008-11-06 10:58:23  rotor
 *  * fixed includes for iostream and friends
 *  * updated for new release (1.0.2)
 *
 * Revision 1.1  2003/05/30 16:43:11  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.5  1996/05/29  16:13:56  rkwan
 * Release 2.5
 *
 * Revision 1.1  1995/12/11  15:29:43  rkwan
 * Initial revision
 *
 *========================================================================*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

/*****************************************************************************
 * Constants
 *****************************************************************************/

// GAMMA_2PI in Hz/T for H-1 nuclei
#define GAMMA_2PI_HZ  42.573E6
#define WATER_FAT_SHIFT_HZ 225

#ifndef TRUE
# define TRUE 1
# define FALSE 0
#endif

#ifndef NULL
# define NULL 0
#endif

#ifdef DEBUG
#include <assert.h>
#include <iostream.h>
#endif

/*****************************************************************************
 * Macros
 *****************************************************************************/

#define DEG_TO_RAD(x) ((double)(M_PI/180.0)*x)
#define RAD_TO_DEG(x) ((double)(180.0/M_PI)*x)

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

/*****************************************************************************
 * Types
 *****************************************************************************/

typedef double Radians;  // angle in Radians
typedef double Degrees;  // angle in Degrees
typedef double Time_ms;  // time in milliseconds
typedef double Time_us;  // time in microseconds
typedef double Hertz;    // frequency in Hz 

typedef int Boolean;  
typedef unsigned char Byte;

#endif

