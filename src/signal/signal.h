#ifndef __SIGNAL_H
#define __SIGNAL_H

/*****************************************************************************
 * 
 * SIGNAL.H
 *
 * MRI Simulator Signal classes
 *
 * R. Kwan
 * March 15, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

#include <math.h>
#include "../mrisim/mrisim.h"

// Magnetization models
//#include "../signal/magvec.h"
//#include "../signal/magsys.h"
#include "../signal/spin_model.h"
#include "../signal/quick_model.h"
#include "../signal/vector_model.h"
#include "../signal/isochromat_model.h"
#include "../signal/fast_iso_model.h"

// Pulse sequences
#include "../signal/pulseseq.h"
#include "../signal/quickseq.h"
#include "../signal/se.h"
#include "../signal/ir.h"
#include "../signal/ce_fast.h"
#include "../signal/fisp.h"
#include "../signal/flash.h"
#include "../signal/spoiled_flash.h"

// Pulse Sequence components
#include "../signal/customseq.h"
#include "../signal/rf_pulse.h"
#include "../signal/sample.h"
#include "../signal/repeat.h"
#include "../signal/spoiler.h"

#include "../signal/tissue.h"

#endif

