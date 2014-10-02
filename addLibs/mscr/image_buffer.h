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
#ifndef __IMAGEBUFFER
#define __IMAGEBUFFER

#define SQR(X) ((X)*(X))

/*   Datatypes for image buffers    */ 

#ifdef DOUBLE_FLOATS
typedef double fpnum;
#else
typedef float fpnum;
#endif

/* Double buffer */
typedef struct buffer {
  int rows;
  int cols;
  int ndim;
  fpnum *data;
} buffer;

/* Integer buffer */
typedef struct ibuffer {
  int rows;
  int cols;
  int ndim;
  int *data;
} ibuffer;

/* Byte buffer */
typedef struct bbuffer {
  int rows;
  int cols;
  int ndim;
  unsigned char *data;
} bbuffer;

/*   Methods for image buffers   */

buffer *buffer_new(int r,int c,int nd);
buffer *buffer_new0(fpnum *darr,int r,int c,int nd);
void    buffer_pdims(buffer *bf);
void    buffer_free(buffer *bf);
buffer *buffer_clone(buffer *bf);

ibuffer *ibuffer_new(int r,int c,int nd);
ibuffer *ibuffer_new0(int *darr,int r,int c,int nd);
void     ibuffer_pdims(ibuffer *bf);
void     ibuffer_free(ibuffer *bf);

bbuffer *bbuffer_new(int r,int c,int nd);
bbuffer *bbuffer_new0(unsigned char *darr,int r,int c,int nd);
void     bbuffer_pdims(bbuffer *bf);
void     bbuffer_free(bbuffer *bf);

#endif /* __IMAGEBUFFER */
