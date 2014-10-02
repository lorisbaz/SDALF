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

#ifndef _MSR_UTIL
#define _MSR_UTIL

#include "image_buffer.h"
#include <limits.h>

#define PI 3.141592653589793143328792211

/* Edge element */
#ifdef INTEGER_EDGES
typedef int edgeval;
#define EDGE_SCALE (INT_MAX/3)
#define EDGE_OFFSET 0.5
#define ebuffer ibuffer
#define ebuffer_encapsulate ibuffer_encapsulate
#define ebuffer_new ibuffer_new
#define ebuffer_free ibuffer_free
typedef struct edge {
  edgeval  dist;
  int  xpos;
  int  ypos;
  int  dirfl;
} edge;
#else
typedef fpnum edgeval;
#define EDGE_SCALE 1
#define EDGE_OFFSET 0
#define ebuffer buffer
#define ebuffer_encapsulate buffer_encapsulate
#define ebuffer_new buffer_new
#define ebuffer_free buffer_free
typedef struct edge {
  edgeval  dist;
  fpnum  xpos;
  fpnum  ypos;
  fpnum  dirfl;
} edge;
#endif

/* bbox element */
typedef struct boundingbox {
  int xmin;
  int xmax;
  int ymin;
  int ymax;
} boundingbox;

/* mvec element */
typedef struct momentvector {
  fpnum area;
  fpnum mx;
  fpnum my;
  fpnum mx2;
  fpnum mxy;
  fpnum my2;
} momentvector;

/* arate element */
typedef struct arearate {
  fpnum arate;
  fpnum amin;
  fpnum d0;
  fpnum dn;
  fpnum t0;
  fpnum tn;
} arearate;

/* dalist element */
typedef struct da_tuple {
  fpnum dist;
  fpnum area;
  int ts;
  int oind;
} da_tuple;

/* node element (doubly linked list) */
typedef struct _node {
  struct _node *prev;
  struct _node *next;
} node;

int edgelist_partition(edge *elist,int nofedges,edgeval pivot);
edgeval edgelist_find_max(ebuffer *bf_elist);
void edgelist_sort(ebuffer *bf_elist);
void filter_print(buffer *bf_filter,char *fname);
buffer *binfilt1d(int order);
buffer *binfilt_der(int order);
buffer *gaussfilt1d(int order);
void blur_buffer(buffer *bf_image0,int order);
void image_colourspace(buffer *bf_img0,buffer *bf_img1);
edgeval image_to_edgelist_grad(buffer *bf_image,ebuffer *bf_elist,int bfz);
edgeval image_to_edgelist_blur(buffer *bf_image,ebuffer *bf_elist,int bfz,int verbosefl);
edgeval image_to_edgelist_blur_n8(buffer *bf_image,ebuffer *bf_elist,int bfz,int verbosefl);
edgeval image_to_edgelist_blur_n8_norm(buffer *bf_image,ebuffer *bf_elist,int bfz,int verbosefl);
edgeval image_to_edgelist_blur_norm(buffer *bf_image,ebuffer *bf_elist,int bfz,int verbosefl);
edgeval image_to_edgelist_blur8(buffer *bf_image,buffer *bf_elist,int bfz);
fpnum image_to_edgelist0(buffer *bf_image,buffer *bf_elist);
void update_edgelist(buffer *bf_image,edge *elist,int nedge,
		     ibuffer *bf_limg,buffer *bf_pvec,buffer *bf_mvec);
void bloblist_defrag(buffer *bf_mvec,buffer *bf_pvec,ibuffer *bf_bbox,
		     buffer *bf_dalist,ibuffer *bf_limg,int *pnbl,int *plno);
void edgelist_to_bloblist(buffer **bfp_mvec,buffer **bfp_pvec,
			  buffer **bfp_arate,buffer *bf_image,
			  ebuffer *bf_elist,ebuffer *bf_thres,int amin,
			  double size_inc,fpnum res,int verbosefl);
void edgelist_to_bloblist_masked(buffer **bfp_mvec,buffer **bfp_pvec,
			  buffer **bfp_arate,buffer *bf_image, buffer *bf_mask,
			  ebuffer *bf_elist,ebuffer *bf_thres,int amin,
			  double size_inc,fpnum res,int verbosefl);
void edgelist_to_bloblist0(buffer **bfp_mvec,buffer **bfp_pvec,
			  buffer **bfp_arate,buffer *bf_image,buffer *bf_elist,
			  buffer *bf_thres,int amin,double size_inc,
			  double defrag_frac,int defrag_min,
			  fpnum res);
void edgelist_to_bloblist00(buffer **bfp_mvec,buffer **bfp_pvec,
			  buffer **bfp_arate,buffer *bf_image,buffer *bf_elist,
			  buffer *bf_thres,int amin,double size_inc,
			  double defrag_frac,int defrag_min,
			  fpnum res);
void edgelist_to_cdf(buffer *bf_elist,buffer *bf_ecdf,edgeval d_max);
void evolution_thresholds0(buffer *bf_ecdf,buffer *bf_thres);
double evolution_thresholds(ebuffer *bf_elist,ebuffer *bf_thres);
double evolution_thresholds2(ebuffer *bf_elist,ebuffer *bf_thres,int order);
void center_moments(buffer *bf_mvec,buffer *bf_pvec);
int bloblist_mark_invalid(buffer *bf_mvec,int amin,buffer *bf_arate,
			  fpnum min_margin);
int bloblist_shape_invalid(buffer *bf_mvec);
void bloblist_compact(buffer *bf_mvec,buffer *bf_mvec2,
		      buffer *bf_pvec,buffer *bf_pvec2,
		      buffer *bf_arate,buffer *bf_arate2);
void bloblist_compact2(buffer *bf_mvec,buffer *bf_mvec2,
		      buffer *bf_pvec,buffer *bf_pvec2);
#endif /* _MSR_UTIL */
