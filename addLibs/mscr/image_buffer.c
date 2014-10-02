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

#include <stdio.h>
#include <stdlib.h>
#include "image_buffer.h"
#include <errno.h>

buffer *buffer_new(int r,int c,int nd) {
buffer *bf;
 fpnum *darr;
 errno=0;
 bf=(buffer *)calloc(1,sizeof(buffer));
 darr=(fpnum *)calloc(r*c*nd,sizeof(fpnum));
 if(errno==ENOMEM) {
   printf("buffer_new: Failed to allocate %dx%dx%dx%ld bytes.\n",
	  r,c,nd,(long)sizeof(fpnum));
 }
 bf->rows=r;
 bf->cols=c;
 bf->ndim=nd;
 bf->data=darr;
 return(bf);
}

/* Version which does not allocate data buffer */
buffer *buffer_new0(fpnum *darr,int r,int c,int nd) {
buffer *bf;
 bf=(buffer *)calloc(1,sizeof(buffer));
 bf->rows=r;
 bf->cols=c;
 bf->ndim=nd;
 bf->data=darr;
 return(bf);
}

/* Make an identical copy of a buffer */

buffer *buffer_clone(buffer *bf) {
int x,y,k;
int rows,cols,ndim,cind;
buffer *bf_new;
fpnum *new_data,*data;
 rows=bf->rows; 
 cols=bf->cols; 
 ndim=bf->ndim;
 data=bf->data;
 /* Allocate new buffer */
 bf_new=buffer_new(rows,cols,ndim);
 new_data=bf_new->data;
 /* Copy buffer, element by element */
 for(x=0;x<cols;x++) {
   for(y=0;y<rows;y++) {
     for(k=0;k<ndim;k++) {
       cind=x*rows+y+k*rows*cols;
       new_data[cind]=data[cind];
     }			  
   }
 }
 return(bf_new);
}

void buffer_pdims(buffer *bf) {
 printf("Buffer: %dx%dx%d\n",bf->rows,bf->cols,bf->ndim);
}

void buffer_free(buffer *bf) {
 free(bf->data);
 free(bf);
}

/*
** Methods for ibuffer
*/

ibuffer *ibuffer_new(int r,int c,int nd) {
ibuffer *bf;
int *darr;
long long nofints;
 errno=0;
 bf=(ibuffer *)calloc(1,sizeof(ibuffer));
 nofints=(long long)r*c*nd;
 darr=(int *)calloc(nofints,sizeof(int));
 if(errno==ENOMEM) {
   printf("ibuffer_new: Failed to allocate %lldx%ld=%dx%dx%dx%ld bytes.\n",
	  nofints,(long)sizeof(int),r,c,nd,(long)sizeof(int));
 }
 bf->rows=r;
 bf->cols=c;
 bf->ndim=nd;
 bf->data=darr;
 return(bf);
}

/* Version which does not allocate data buffer */
ibuffer *ibuffer_new0(int *darr,int r,int c,
		      int nd) {
ibuffer *bf;
 bf=(ibuffer *)calloc(1,sizeof(ibuffer));
 bf->rows=r;
 bf->cols=c;
 bf->ndim=nd;
 bf->data=darr;
 return(bf);
}

void ibuffer_pdims(ibuffer *bf) {
 printf("iBuffer: %dx%dx%d\n",bf->rows,bf->cols,bf->ndim);
}

void ibuffer_free(ibuffer *bf) {
 free(bf->data);
 free(bf);
}

/*
** Methods for bbuffer
*/

bbuffer *bbuffer_new(int r,int c,int nd) {
bbuffer *bf;
unsigned char *darr;
long long nofbytes;
 errno=0;
 bf=(bbuffer *)calloc(1,sizeof(bbuffer));
 nofbytes=(long long)r*c*nd;
 darr=(unsigned char *)calloc(nofbytes,sizeof(unsigned char));
 if(errno==ENOMEM) {
   printf("bbuffer_new: Failed to allocate %lldx%ld=%dx%dx%dx%ld bytes.\n",
	  nofbytes,(long)sizeof(unsigned char),r,c,nd,(long)sizeof(unsigned char));
 }
 bf->rows=r;
 bf->cols=c;
 bf->ndim=nd;
 bf->data=darr;
 return(bf);
}

/* Version which does not allocate data buffer */
bbuffer *bbuffer_new0(unsigned char *darr,int r,int c,
		      int nd) {
bbuffer *bf;
 bf=(bbuffer *)calloc(1,sizeof(bbuffer));
 bf->rows=r;
 bf->cols=c;
 bf->ndim=nd;
 bf->data=darr;
 return(bf);
}

void bbuffer_pdims(bbuffer *bf) {
 printf("bBuffer: %dx%dx%d\n",bf->rows,bf->cols,bf->ndim);
}

void bbuffer_free(bbuffer *bf) {
 free(bf->data);
 free(bf);
}

