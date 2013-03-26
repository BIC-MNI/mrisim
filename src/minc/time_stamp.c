/* ----------------------------- MNI Header -----------------------------------
@NAME       : time_stamp.c
@DESCRIPTION: File containing routine to create a time stamp string.
@METHOD     : 
@CREATED    : February 1, 1993 (Peter Neelin)
@MODIFIED   : $Log: time_stamp.c,v $
@MODIFIED   : Revision 1.1  2003-05-30 16:43:10  bert
@MODIFIED   : Initial checkin, mrisim 3.1 from Remi Kwan's home directory
@MODIFIED   :
 * Revision 2.5  1996/05/29  16:30:04  rkwan
 * Release 2.5
 *
 * Revision 1.1  1995/10/13  14:19:23  rkwan
 * Initial revision
 *
 * Revision 1.1  1995/07/26  15:04:49  rkwan
 * Initial revision
 *
 * Revision 1.3  93/08/04  13:03:56  neelin
 * Added RCS $Log: time_stamp.c,v $
 * Added RCS Revision 1.1  2003-05-30 16:43:10  bert
 * Added RCS Initial checkin, mrisim 3.1 from Remi Kwan's home directory
 * Added RCS
 * Revision 2.5  1996/05/29  16:30:04  rkwan
 * Release 2.5
 *
 * Revision 1.1  1995/10/13  14:19:23  rkwan
 * Initial revision
 *
 * Revision 1.1  1995/07/26  15:04:49  rkwan
 * Initial revision
 * to keep track of modifications in source.
 * 
@COPYRIGHT  :
              Copyright 1993 Peter Neelin, McConnell Brain Imaging Centre, 
              Montreal Neurological Institute, McGill University.
              Permission to use, copy, modify, and distribute this
              software and its documentation for any purpose and without
              fee is hereby granted, provided that the above copyright
              notice appear in all copies.  The author and McGill University
              make no representations about the suitability of this
              software for any purpose.  It is provided "as is" without
              express or implied warranty.
---------------------------------------------------------------------------- */
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "time_stamp.h"

#ifndef lint
static char rcsid[]="$Header: /static-cvsroot/simulation/mrisim/src/minc/time_stamp.c,v 1.1 2003-05-30 16:43:10 bert Exp $";
#endif

/* ----------------------------- MNI Header -----------------------------------
@NAME       : time_stamp
@INPUT      : argc - number of arguments
              argv - list of arguments
@OUTPUT     : 
@RETURNS    : pointer to string containing time stamp.
@DESCRIPTION: Function to produce a time stamp string for a program.
              Returns a string of the form "date > command". The command
              is simply the concatenation of argv elements.
@METHOD     : 
@GLOBALS    : 
@CALLS      : 
@CREATED    : February 1, 1993 (Peter Neelin)
@MODIFIED   : 
---------------------------------------------------------------------------- */
char *time_stamp(int argc, char *argv[])
{
   char *str, *the_time;
   int length, i, last;
   static char separator[]={">>>"};
   time_t timer;

   /* Get the time, overwriting newline */
   timer = time(NULL);
   the_time = ctime(&timer);

   /* Get the total length of the string and allocate space */
   length=strlen(the_time) + strlen(separator) + 2;
   for(i=0; i<argc; i++) {
      length += strlen(argv[i]) + 1;
   }
   str = malloc(length);

   /* Copy the time and separator */
   (void) strcpy(str, the_time);
   str[strlen(str)-1]='\0';
   (void) strcat(str, separator);

   /* Copy the program name and arguments */
   for (i=0; i<argc; i++) {
      last = strlen(str);
      str[last]=' ';
      str[last+1]='\0';
      (void) strcat(str, argv[i]);
   }

   /* Add a terminating newline */
   last = strlen(str);
   str[last]='\n';
   str[last+1]='\0';


   return str;
}
