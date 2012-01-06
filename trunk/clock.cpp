//  file CLOCK.C

/*
  This file customizes function 'timeNow' to your platform.
  Function 'timeNow' is used for measuring elapsed or CPU time.
*/

#include "clock.h"


#ifdef BORLAND
#include <time.h>
   float timeNow () {return clock()/CLOCKS_PER_SEC ;}
#endif

#ifdef SILICON_GRAPHICS
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
   float timeNow ()
   {
      float       t,sum ;
      struct tms  buffer ;

      times(&buffer) ;
      sum = (buffer.tms_utime) + (buffer.tms_stime) ;
      t = sum/HZ ;
      // printf ("Total user and CPU time %f\n",t) ;
      return  t ;
   }
#endif

#ifdef SUN_STATION
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>
   float timeNow ()
   {
      float       t,sum ;
      struct tms  buffer ;

      times(&buffer) ;
      sum = (buffer.tms_utime) + (buffer.tms_stime) ;
      t = sum/HZ ;
      // printf ("Total user and CPU time %f\n",t) ;
      return  t ;
   }
#endif

#ifdef OTHER_PLATFORM
   float timeNow () {return 0. ;}
#endif

