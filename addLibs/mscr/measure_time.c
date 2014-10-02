/*
**  Source file for mscr. (c) 2007 Per-Erik Forssen
**
**  This program is free software; you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation; either version 2 of the License, or
**  (at your option) any later version.
**  
**  See the file COPYING for details.
**
*/

#ifndef WIN32
#include <sys/time.h>
#endif

#include <stdlib.h>

#include "measure_time.h"

/*
 *  Time, accurate up to microseconds
 */
double accurate_clock() {
#ifndef WIN32
  struct timeval tp;
  gettimeofday(&tp,NULL);

  return tp.tv_sec + tp.tv_usec/1e6;
#else
  return 0
#endif
}

timewin *timewin_new(int tsize) {
#ifndef WIN32
  timewin *tw;
  tw=calloc(1,sizeof(timewin));
  tw->nticks=tsize;
  tw->ticks=calloc(tw->nticks,sizeof(double));
  return tw;
#else
  return NULL
#endif
}

void timewin_free(timewin *tw) {
#ifndef WIN32
  free(tw->ticks);
  free(tw);
#endif
}

timewin *timewin_addtime(timewin *tw) {
#ifndef WIN32
  int k;
  /* Shift old times */
  for(k=1;k<tw->nticks;k++) {
    tw->ticks[k-1]=tw->ticks[k];
  }
  tw->ticks[tw->nticks-1]=accurate_clock();  /* Add new time */
  return tw;
#else
  return 0;
#endif
}

/*
**  Compute rate (typically framerate) from a window of measurements
*/
double timewin_rate(timewin *tw) {
#ifndef WIN32
  return (double)(tw->nticks-1.0)/(tw->ticks[tw->nticks-1]-tw->ticks[0]);
#else
  return 0;
#endif
}
