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
#include <string.h>
#endif
#include "mex.h"
#include "mexutil.h"
#include "image_buffer.h"
#include "visualise.h"

/*
**  mexFunction is the interface function to MATLAB.
**
** plhs - pointer to left-hand OUTPUT mxArrays
** nlhs - number of left-hand OUTPUT arrays
** prhs - pointer to right-hand INPUT mxArrays
** nrhs - number of right-hand INPUT arrays
**
*/
void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
fpnum bkgr_green[] = {0.0,1.0,0.0}; /* Background colour */
buffer *bf_mvec,*bf_pvec,*bf_pvecn,*bf_rimg,*bf_blobimg;
buffer *bf_bkgr;
ibuffer *bf_limg;
const int *arg_dims;
int arg_ndims,grxfl,rows,cols;
int out_dims[]={0,0,0};
int mdims[]={0,0,0};

/* Argument count check */
 if((nrhs <4)||(nrhs>5)) {
   if(nrhs !=0) mexPrintf("Error: 4-5 input arguments expected.\n");
   goto mexErrExit;
 }
    
/* Parse arg 1 (mvec) */
   mdims[0]=6;mdims[1]=-1;
   if(!type_check("mvec",prhs[0],mxDOUBLE_CLASS,2,mdims))
     goto mexErrExit;
   bf_mvec=buffer_encapsulate(prhs[0]);

/* Parse arg 2 (pvec) */
   mdims[0]=-1;mdims[1]=bf_mvec->cols;
   if(!type_check("pvec",prhs[1],mxDOUBLE_CLASS,2,mdims))
     goto mexErrExit;
   bf_pvec=buffer_encapsulate(prhs[1]);

/* Parse arg 3 (rows) */
   mdims[0]=1;mdims[1]=1;
   if(!type_check("rows",prhs[2],mxDOUBLE_CLASS,2,mdims))
     goto mexErrExit;
   rows=mxGetScalar(prhs[2]);

/* Parse arg 4 (cols) */
   mdims[0]=1;mdims[1]=1;
   if(!type_check("cols",prhs[3],mxDOUBLE_CLASS,2,mdims))
     goto mexErrExit;
   cols=mxGetScalar(prhs[3]);

   /* Set image size */
   out_dims[0]=rows;
   out_dims[1]=cols;

   /* Copy number of colour bands */
   out_dims[2]=bf_pvec->rows;

   /* Parse arg 4 (bkgr) */
   if(nrhs==5) {
     mdims[0]=out_dims[2];mdims[1]=1;
     if(!type_check("bkgr",prhs[4],mxDOUBLE_CLASS,2,mdims))
       goto mexErrExit;
     bf_bkgr=buffer_encapsulate(prhs[4]);
 } else {
   bf_bkgr=buffer_new0(bkgr_green,3,1,1);  /* Use default */
 }

 /* Extend grey to RGB? */
 grxfl=0;
 if(out_dims[2]==1) {
   grxfl=1;
   pvec_to_rgb(bf_pvec,&bf_pvecn);
   out_dims[2]=3;
 }
/* Call visualisation functions */
 if(nlhs>0) {
   /* Create an empty green image */
   plhs[0] = mxCreateNumericArray(3,out_dims,mxDOUBLE_CLASS,mxREAL);
   bf_blobimg = buffer_encapsulate(plhs[0]);
   buffer_paint(bf_blobimg,bf_bkgr);
   
   /* Set area to approximating ellipse area */
   mvec_set_area(bf_mvec);  /* Blobs are drawn in descending area order */

   /* Visualise blobs in the green image */
   if(grxfl)
     draw_ellipses(bf_blobimg,bf_mvec,bf_pvecn);
   else
     draw_ellipses(bf_blobimg,bf_mvec,bf_pvec);
   free(bf_blobimg); /* Only release struct */
 }

/* Free memory */
 free(bf_mvec);  /* Only release struct */
 free(bf_bkgr);  /* Only release struct */
 free(bf_pvec);  /* Only release struct */
 if(grxfl) buffer_free(bf_pvecn);
     
/* Jump to normal exit */
 goto mexExit;

/* Exit with error message */
mexErrExit:
 mexPrintf("Usage: <Bimg>=draw_blobs(<mvec>,<pvec>,<rows>,<cols>[,<bkgr>])\n");
mexExit: {}
}
