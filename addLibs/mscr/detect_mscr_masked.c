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
**
**  Added mask foreground/background 
**  Michela Farenzena, Sept 2009
 */

#include <math.h>
#include <string.h>
#include "mex.h"
#include "image_buffer.h"
#include "mexutil.h"
#include "msr_util.h"

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
int msize[]={0,0,0};
int arg_ndims;
int out_dims[]={0,0,0};
buffer *bf_image,*bf_pvec,*bf_pvec2, *bf_mask;
ebuffer *bf_elist,*bf_thres=NULL;
buffer *bf_mvec,*bf_mvec2,*bf_arate,*bf_arate2;
buffer *bf_ecdf;
int validcnt,nofedges,rows,cols,ndim;
edgeval d_max;
double d_mean;
mxArray *mx_pvec,*mx_arate,*mx_elist;
int tslist_flag=1;
int min_size=60;           /* Default */
double ainc=1.05;          /* Default */
double min_marginf=0.0015; /* Default */
edgeval min_margin;        /* Value after conversion */
fpnum res=1e-4;            /* Default */
int Ns=10;                 /* Default */
int timesteps=200;         /* Default */
int filter_size=3;         /* Default */
int n8flag=0;              /* Default */
int normfl=1;              /* Default */
int blurfl=0;              /* Default */
int verbosefl=1;           /* Default */

   /* Argument count check */
  if((nrhs ==0)||(nrhs > 3)) {
    mexPrintf("Usage: [<mvec>,<pvec>,<arate>,<elist2>]=detect_mscr(<image>,<mask>[,<pars>])\n");
    
    if((nrhs > 3))
      mexErrMsgTxt("Error: 2 or 3 input arguments required.\n");
    return;
  }

  /* Parse arg 1 (RGB or greyscale image) */
  msize[0]=-1;msize[1]=-1;msize[2]=-1;
  arg_ndims = mxGetNumberOfDimensions(prhs[0]);
  if((arg_ndims!=2)&&(arg_ndims!=3)) {
    mexPrintf("Error: <image> should be 2D or 3D.\n");
    return;
  }
  if(!type_check("image",prhs[0],mxDOUBLE_CLASS,arg_ndims,msize))
    return;
  
  bf_image=buffer_encapsulate(prhs[0]);
  rows = bf_image->rows;
  cols = bf_image->cols;
  ndim = bf_image->ndim;
  if((ndim!=1)&&(ndim!=3)) {
    mexPrintf("Size mismatch: <Im> should be MxNxD where D=1 or 3.\n");
    free(bf_image);
    return;
  }
  
  

 

   /* Parse arg 2 (mask image) */ 
  bf_mask=buffer_encapsulate(prhs[1]);
  rows = bf_mask->rows;
  cols = bf_mask->cols;
  ndim = bf_mask->ndim;
  if((ndim!=1)&&(rows!=bf_image->rows)&&cols!=bf_image->cols) {
    mexPrintf("Size mismatch: <mask> should be of the same dimention as <image>.\n");
    free(bf_image);
    return;
  }
  
    
  /* Parse arg 3 (pars struct) */
  if(nrhs>2) {
    msize[0]=1;msize[1]=1;
    if(!type_check("pars",prhs[2],mxSTRUCT_CLASS,2,msize)) {
      free(bf_image);
      return;
    }

    min_marginf = get_field_value(prhs[2],"min_margin",min_marginf);
    min_margin  = min_marginf*EDGE_SCALE+EDGE_OFFSET;
    timesteps   = get_field_value(prhs[2],"timesteps",timesteps);
    min_size    = get_field_value(prhs[2],"min_size",min_size);
    ainc        = get_field_value(prhs[2],"ainc",ainc);
    filter_size = get_field_value(prhs[2],"filter_size",filter_size);
    n8flag      = get_field_value(prhs[2],"n8flag",n8flag);
    normfl      = get_field_value(prhs[2],"normfl",normfl);
    blurfl      = get_field_value(prhs[2],"blurfl",blurfl);
    verbosefl   = get_field_value(prhs[2],"verbosefl",verbosefl);
  }

  if(verbosefl) {
    /* Display parameter settings */
    mexPrintf("\nDescriptor parameter settings:\n");
    mexPrintf(" min_margin: %g\n",min_marginf);
    mexPrintf("  timesteps: %d\n",timesteps);
    mexPrintf("   min_size: %d\n",min_size);
    mexPrintf("       ainc: %g\n",ainc);
    mexPrintf("filter_size: %d\n",filter_size);
    mexPrintf("     n8flag: %d\n",n8flag);
    mexPrintf("     normfl: %d\n",normfl);
    mexPrintf("     blurfl: %d\n",blurfl);
    mexPrintf("  verbosefl: %d\n",verbosefl);
  }
  
 
  
  
 
 /* Set edge list size */
 if(n8flag) 
   nofedges=4*rows*cols-3*rows-3*cols+2;
 else
   nofedges=2*rows*cols-rows-cols;

 if(nlhs>0) {
   /* Allocate out elist */
   out_dims[0]=4;
   out_dims[1]=nofedges;
#ifdef INTEGER_EDGES
   mx_elist = mxCreateNumericArray(2,out_dims,mxINT32_CLASS,mxREAL);
#else
   mx_elist = mxCreateNumericArray(2,out_dims,mxDOUBLE_CLASS,mxREAL);
#endif
   bf_elist = ebuffer_encapsulate(mx_elist);   

   /* Call computation function */

   if(blurfl) {
     /* Blur input image instead */
     blur_buffer(bf_image,filter_size);
     filter_size=1;
   }
   if(filter_size%2) {
     if(n8flag) {
       if(normfl)
	 d_max=image_to_edgelist_blur_n8_norm(bf_image,bf_elist,filter_size,verbosefl);
       else
	 d_max=image_to_edgelist_blur_n8(bf_image,bf_elist,filter_size,verbosefl);
     } else {
       if(normfl)
	 d_max=image_to_edgelist_blur_norm(bf_image,bf_elist,filter_size,verbosefl);
       else
	 d_max=image_to_edgelist_blur(bf_image,bf_elist,filter_size,verbosefl);
     }
   } else {
     mexErrMsgTxt("bfz should be odd.\n");
/*      d_max=image_to_edgelist_grad(bf_image,bf_elist,filter_size); */
   }


   if(tslist_flag) {
     /* Call cdf sampling */
/*      bf_ecdf=buffer_new(2,Ns,1); */
/*      edgelist_to_cdf(bf_elist,bf_ecdf,d_max); */

     /* Call thresholds interpolation */
     bf_thres=ebuffer_new(1,timesteps,1);
     if(verbosefl) mexPrintf("order=%d\n",ndim);
     d_mean=evolution_thresholds2(bf_elist,bf_thres,ndim);
     if(verbosefl) mexPrintf("d_mean=%g\n",d_mean);
     /* Try linear dependence on mean edge strength */
/*      min_margin=d_mean*min_margin_sc; */
   }

   /* This is the main function */
   edgelist_to_bloblist_masked(&bf_mvec,&bf_pvec,&bf_arate,bf_image,bf_mask,bf_elist,bf_thres,min_size,ainc,res,verbosefl);
   center_moments(bf_mvec,bf_pvec);
   validcnt=bloblist_mark_invalid(bf_mvec,min_size,bf_arate,(fpnum)min_margin);
   validcnt=bloblist_shape_invalid(bf_mvec);
   if(verbosefl) mexPrintf("validcnt=%d\n",validcnt);

   if(tslist_flag) {
     ebuffer_free(bf_thres);
/*      buffer_free(bf_ecdf); */
   }

   /* Allocate out mvec */
   out_dims[0]=6;
   out_dims[1]=validcnt;
   plhs[0] = mxCreateNumericArray(2,out_dims,mxDOUBLE_CLASS,mxREAL);
   bf_mvec2 = buffer_encapsulate(plhs[0]);

   /* Allocate out pvec */
   out_dims[0]=3;
   out_dims[1]=validcnt;
   mx_pvec = mxCreateNumericArray(2,out_dims,mxDOUBLE_CLASS,mxREAL);
   bf_pvec2 = buffer_encapsulate(mx_pvec);

   /* Allocate out arates */
   out_dims[0]=bf_arate->rows;
   out_dims[1]=validcnt;
   mx_arate = mxCreateNumericArray(2,out_dims,mxDOUBLE_CLASS,mxREAL);
   bf_arate2 = buffer_encapsulate(mx_arate);

   bloblist_compact(bf_mvec,bf_mvec2,bf_pvec,bf_pvec2,bf_arate,bf_arate2);
   free(bf_mvec2);    /* Only release struct */
   free(bf_pvec2);    /* Only release struct */
   free(bf_arate2);   /* Only release struct */
   free(bf_elist);    /* Only release struct */

   buffer_free(bf_mvec);   /* Release non-compacted mvec */
   buffer_free(bf_pvec);   /* Release non-compacted pvec */
   buffer_free(bf_arate);  /* Release non-compacted arate */
 }

 if(nlhs>1) plhs[1]=mx_pvec; 

 if(nlhs>2) plhs[2]=mx_arate; 

 if(nlhs>3) plhs[3]=mx_elist; 

 /* Free memory */
 free(bf_image);  /* Only release struct */
}
