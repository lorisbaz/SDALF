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
#ifndef __VISUALISE
#define __VISUALISE

#include "image_buffer.h"
#define PI 3.141592653589793143328792211

/*   Visualisation methods  */

void buffer_paint(buffer *bf,buffer *bf_pvec);
void bbuffer_paint(bbuffer *bf,bbuffer *bf_pvec);
void eigendec(fpnum *I,fpnum *D,fpnum *E);
void draw_ellipses(buffer *bf_img,buffer *bf_mvec,buffer *bf_pvec);
void bdraw_ellipses(bbuffer *bf_img,buffer *bf_mvec,bbuffer *bf_pvec);
void draw_regions(buffer *bf_img,ibuffer *bf_labelim,buffer *bf_pvec);
void pvec_to_rgb(buffer *bf_pvec,buffer **bl_pvecn);
void pvec_colourspace(buffer *bf_pvec);
void pvec_to_uint8(buffer *bf_pvec,bbuffer **bl_pvec2);
void mvec_set_area(buffer *bf_mvec);
void inside_mask_ratios(ibuffer *bf_mask,buffer *bf_mvec,buffer *bf_ratios);

#endif /* __VISUALISE */
