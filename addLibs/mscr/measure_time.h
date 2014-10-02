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

#ifndef _MEASURE_TIME
#define _MEASURE_TIME

typedef struct timewin {
  int nticks;
  double *ticks;
} timewin;

double accurate_clock();

timewin *timewin_new(int tsize);
void timewin_free(timewin *tw);
timewin *timewin_addtime(timewin *tw);
double timewin_rate(timewin *tw);

#endif /* _MEASURE_TIME */
