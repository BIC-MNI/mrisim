#ifndef __CUSTOMSEQ_H
#define __CUSTOMSEQ_H

/*****************************************************************************
 * 
 * CUSTOMSEQ.H
 *
 * Custom_Sequence Class
 *
 * R. Kwan
 * August 2, 1995
 *
 * (C) Copyright 1995 by R.Kwan
 *
 *****************************************************************************/

/* $Header: /private-cvsroot/simulation/mrisim/src/signal/customseq.h,v 1.1 2003-05-30 16:43:12 bert Exp $
 * $Log: customseq.h,v $
 * Revision 1.1  2003-05-30 16:43:12  bert
 * Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 *
 * Revision 2.3  1995/12/11  14:22:56  rkwan
 * Updated for fast_iso_model.
 *
 * Revision 2.2  1995/08/08  15:26:57  rkwan
 * Eliminated #include's.
 *
 * Revision 2.1  1995/08/08  14:32:13  rkwan
 * Updated for new spin model and pulse sequence.
 *
 */

#include <mrisim/mrisim.h>

#include "pulseseq.h"
#include "event.h"
#include "spin_model.h"

/*****************************************************************************
 * Custom_Sequence Class
 *****************************************************************************/

class Custom_Sequence : public Pulse_Sequence {
   public:
      Custom_Sequence();                        // constructors
      Custom_Sequence(const Pulse_Sequence& p);
      Custom_Sequence(const Custom_Sequence& p);

      virtual ~Custom_Sequence();               // destructor

      // Pulse Sequence programming
      void add_event(Event *event);
      void add_event(Time_ms t, Event *event);
      int  delete_event(int event_id);
      int  get_next_event_id(void);

      // Pulse Sequence acquisition
      Vector_3D& initialize_sequence(Spin_Model& model);
      Vector_3D& apply_to_time(Spin_Model& model, Time_ms t);
      Vector_3D& apply_next_event(Spin_Model& model);
      Vector_3D& apply_to_end_of_repetition(Spin_Model& model);
      Vector_3D& apply_one_repetition(Spin_Model& model);
      Vector_3D& apply_to_steady_state(Spin_Model& model);

      void apply_trace(Spin_Model& model, int trace_length, 
                       Time_ms trace_step, Time_ms time[],
                       Vector_3D m[]);

      void display_sequence_info(ostream& stream);
      void dump_sequence_info(FILE *output);

   protected:

      int     _num_events;
      Event   *_event_list;
      Event   *_next_event;

      void    _delete_list(Event *head);
      Event   *_find_event_id(int event_id);
      Event   *_find_event_after_time(Time_ms t);
};

#endif



