#ifndef __MRISIM_MAIN_H
#define __MRISIM_MAIN_H

//==========================================================================
// MRISIM_MAIN.H
// Include header and function prototypes for command line version
// of mrisim.
// Inherits from:
// Base class to:  
//
// R. Kwan
// (C) Copyright 1995 by R.Kwan
//==========================================================================

/*==========================================================================
 * $Header: /private-cvsroot/simulation/mrisim/src/mrisim/mrisim_main.h,v 1.1 2003-05-30 16:43:11 bert Exp $
 * $Log: mrisim_main.h,v $
 * Revision 1.1  2003-05-30 16:43:11  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 3.1  1996/07/19  15:58:25  rkwan
 * Release 3.1 update.
 *
 * Revision 2.5  1996/05/29  16:15:33  rkwan
 * Release 2.5
 *
 * Revision 1.6  1996/01/17  18:02:11  rkwan
 * Update for tx/rx coil modelling.
 *
 * Revision 1.5  1995/12/22  20:16:52  rkwan
 * Update for percent coil and RF map features.
 *
 * Revision 1.4  1995/12/11  17:50:17  rkwan
 * Check for output file existence added.
 *
 * Revision 1.3  1995/12/11  15:59:02  rkwan
 * RCS fix.
 *
 * Revision 1.2  1995/12/11  15:24:34  rkwan
 * Doc update.
 *
 *========================================================================*/

#include "mriscanner.h"
#include "mrisimargs.h"
#include "paramfile.h"

#include "discrete_phantom.h"
#include "discrete_rf_phantom.h"
#include "fuzzy_phantom.h"
#include "fuzzy_rf_phantom.h"

#include "rf_coil.h"
#include "intrinsic_coil.h"
#include "image_snr_coil.h"
#include "percent_coil.h"

#include <signal/signal.h>
#include <minc/mriminc.h>

void create_models(const mrisimArgs &args, MRI_Scanner &scanner);

void output_images(const mrisimArgs &args, char *stamp, MRI_Scanner &scanner);

Phantom *open_phantom(const mrisimArgs &args, const RF_Coil *rf_coil);

RF_Coil *make_coil(const mrisimArgs &args);

int apply_pulse_sequence(const mrisimArgs &args, MRI_Scanner &scanner);

#endif
