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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "image_buffer.h"
#include "msr_util.h"

/*
** Partitioning function
** Partitions the data according to <pivot>
** and returns number of smaller elements
*/
int edgelist_partition(edge *elist,int nofedges,edgeval pivot) {
int tpos,inspos;
edge etemp;
 inspos=0;
 for(tpos=0;tpos<nofedges;tpos++) {
   if(elist[tpos].dist<pivot) {
     etemp=elist[tpos];
     elist[tpos]=elist[inspos];
     elist[inspos++]=etemp;
   }
 }
 return inspos;
}
/*
 * Sort elements of an edgelist in increasing order
 * Using quicksort. Unpredictable time when many
 * elements are similar.
 */
void elist_sort(edge *elist,int nofedges) {
edgeval p1,p2,p3,pivot;
/* edgeval minval; */
edge etemp;
int lpart,k,l;
/* int minind; */
 if(nofedges<10) {
   /* Sort in place... (insertion sort) */
   for(k=1;k<nofedges;k++) {
     etemp=elist[k];
     l=k-1;
     while((l>=0)&&(elist[l].dist>etemp.dist)) {
       elist[l+1]=elist[l];
       l--;
     }
     if(l<k-1) elist[l+1]=etemp;
   }

   /* Sort in place... (selection sort) */
/*    for(k=0;k<nofedges;k++) { */
/*      minval=elist[k].dist; */
/*      minind=k; */
/*      for(l=k+1;l<nofedges;l++) { */
/*        if(elist[l].dist<minval) { */
/* 	 minval=elist[l].dist; */
/* 	 minind=l; */
/*        } */
/*      } */
/*      /\* Swap *\/ */
/*      if(minind>k) { */
/*        etemp=elist[k]; */
/*        elist[k]=elist[minind]; */
/*        elist[minind]=etemp; */
/*      } */
/*    } */
 } else {
   /* Find pivot as median of 3 */
   p1=elist[0].dist;
   p2=elist[nofedges/2].dist;
   p3=elist[nofedges-1].dist;
/*    p2=elist[1].dist; */
/*    p3=elist[2].dist; */
   if(p1>p2) {
     if(p1>p3) {
       if(p2>p3) pivot=p2; else pivot=p3;
     } else {
       pivot=p1;
     }
   } else {
     if(p2>p3) {
       if(p3>p1) pivot=p3; else pivot=p1;
     } else {
       pivot=p2;
     }
   }
   /* Randomised pivot */
/*    pivot=elist[(rand()*nofedges)/RAND_MAX].dist; */
   /* */
   lpart=edgelist_partition(elist,nofedges,pivot);
/*    mexPrintf("lpart=%d nofedges=%d pivot=%g\n",lpart,nofedges,(double)pivot); */
   if(lpart>0) {
     elist_sort(&elist[lpart],nofedges-lpart);
     elist_sort(elist,lpart);
   } else {
   /* Sort in place... (insertion sort) */
   for(k=1;k<nofedges;k++) {
     etemp=elist[k];
     l=k-1;
     while((l>=0)&&(elist[l].dist>etemp.dist)) {
       elist[l+1]=elist[l];
       l--;
     }
     if(l<k-1) elist[l+1]=etemp;
   }

     /* Sort in place... (selection sort) */
/*      for(k=0;k<nofedges;k++) { */
/*        minval=elist[k].dist; */
/*        minind=k; */
/*        for(l=k+1;l<nofedges;l++) { */
/* 	 if(elist[l].dist<minval) { */
/* 	   minval=elist[l].dist; */
/* 	   minind=l; */
/* 	 } */
/*        } */
/*        /\* Swap *\/ */
/*        if(minind>k) { */
/* 	 etemp=elist[k]; */
/* 	 elist[k]=elist[minind]; */
/* 	 elist[minind]=etemp; */
/*        } */
/*      } */
   }
 }
}
/*
 * Sort elements of an edgelist in increasing order
 */
void edgelist_sort(ebuffer *bf_elist) {
 elist_sort((edge *)bf_elist->data,bf_elist->cols);
}
/* 
 * Loop through an edge list and find max
 */
edgeval edgelist_find_max(ebuffer *bf_elist) {
int nedges,x;
edge *elist;
edgeval d_max;
 nedges=bf_elist->cols; 
 elist = (edge *)bf_elist->data;
 d_max=0.0;
 for(x=0;x<nedges;x++) {
   if(elist[x].dist>d_max) d_max=elist[x].dist;
 }
 return d_max;
}
/*
 * Print filter coefficients (for debugging)
 */
void filter_print(buffer *bf_filter,char *varname) {
fpnum *filter;
int x,ncoeff;
 ncoeff = bf_filter->rows;
 filter = bf_filter->data;
 printf("%s=[",varname);
 for(x=0;x<ncoeff;x++) {
   printf("%g ",filter[x]);
 }
 printf("]\n");
}
/*
 *  Binfilt function
 *
 *  Allocates space and returns a buffer
 *  containing a binomial filter normalised
 *  to sum to 1.
 *
 */
buffer *binfilt1d(int order) {
buffer *bf_w1,*bf_w2;
fpnum *w1,*w2,*wtmp;
fpnum wsum;
int k,l;
  bf_w1=buffer_new(order,1,1);
  bf_w2=buffer_new(order,1,1);
  w1=bf_w1->data;
  w2=bf_w2->data;
  w1[0]=1;
  w2[0]=1;
  for(k=0;k<order-1;k++) {
    for(l=1;l<order;l++) {
      w2[l]=w1[l-1]+w1[l];
    }
    wtmp = w1;
    w1 = w2;
    w2 = wtmp;
  }
  wsum=exp(((fpnum)order-1.0)*log(2));
  for(k=0;k<order;k++) {
    w1[k]=w1[k]/wsum;
  }
  if(order%2) {
    /* Odd -> result is in bf_w1 */
    buffer_free(bf_w2);
    return bf_w1;
  } else {
    /* Even -> result is in bf_w2 */
    buffer_free(bf_w1);
    return bf_w2;
  }
}
/*
 *  Generate a 1D binfilt derivative filter
 *  (should be even)
 */
buffer *binfilt_der(int order) {
buffer *bf_w1;
fpnum *w1,origin,dxn;
int x;
 bf_w1 = binfilt1d(order);   /* Start off from binfilt */
 w1 = bf_w1->data; 
 origin=((fpnum)order-1.0)/2.0;
 dxn=0.0;
 for(x=0;x<order;x++) {
   w1[x]=w1[x]*((fpnum)x-origin);
   dxn += fabs(w1[x]*((fpnum)x-origin));
/*     dx=dx/(abs(x*dx'));         % dx of ramp should be 1 */
 }
 for(x=0;x<order;x++) {
   w1[x]=w1[x]/dxn;
 }
 return bf_w1;
}
/*
 * Gaussian filter function
 *
 *  Allocates space and returns a buffer
 *  containing a Gaussian filter normalised
 *  to sum to 1.
 *
 */
buffer *gaussfilt1d(int order) {
buffer *bf_w;
fpnum *w;
fpnum std,wsum,mu;
int k;
  bf_w=buffer_new(order,1,1);
  w=bf_w->data;
  mu=(fpnum)(order-1)/2.0;
  std=sqrt((fpnum)order/5);   /* Same std as binomial */
  wsum=0;
  for(k=0;k<order;k++) {
    w[k]=exp(-.5*(k-mu)*(k-mu)/std/std);
    wsum+=w[k];
    }
  for(k=0;k<order;k++) w[k] /=wsum;
  return bf_w;
}
/*
 * Destructively blur image buffer image0
 * with binomial or gaussian filter of given order.
 */
void blur_buffer(buffer *bf_image0,int order) {
int rows,cols,ndim,field_offset,x,y,m,f;
buffer *bf_filter,*bf_image1; 
fpnum *image0,*image1,*filter;
int fmin,fmax;
fpnum cw,ws,bs;
int gaussfl=1; 
 rows  = bf_image0->rows;
 cols  = bf_image0->cols;
 ndim  = bf_image0->ndim;
 bf_image1 = buffer_new(rows,cols,ndim);
 image0 = bf_image0->data;
 image1 = bf_image1->data;
 if(gaussfl) 
   bf_filter = gaussfilt1d(order);
 else
   bf_filter = binfilt1d(order);
 filter = bf_filter->data;
 fmin=-floor(order/2);
 fmax=floor(order/2);
 filter = bf_filter->data;
 for(m=0;m<ndim;m++) {
   field_offset=rows*cols*m;
   for(x=0;x<cols;x++) {
     for(y=0;y<rows;y++) {
       ws=0.0;bs=0.0;
       for(f=fmin;f<=fmax;f++) {
	 if((x+f>=0)&&(x+f<cols)) {
	   cw=filter[f-fmin];
	   ws += cw;
	   bs += cw*image0[y+(x+f)*rows+field_offset];
	 }
       }
       image1[y+x*rows+field_offset]=bs/ws;
     }
   }
   for(x=0;x<cols;x++) {
     for(y=0;y<rows;y++) {
       ws=0.0;bs=0.0;
       for(f=fmin;f<=fmax;f++) {
	 if((y+f>=0)&&(y+f<rows)) {
	   cw=filter[f-fmin];
	   ws += cw;
	   bs += cw*image1[y+f+x*rows+field_offset];
	 }
       }
       image0[y+x*rows+field_offset]=bs/ws;
     }
   }
 }
 buffer_free(bf_image1);
 buffer_free(bf_filter);
}
/*
 * Switch colourspace
 */
void image_colourspace(buffer *bf_img0,buffer *bf_img1) {
  fpnum *img0,*img1,cval;
 int rows,cols,ndim,x,y;
 rows=bf_img0->rows;
 cols=bf_img0->cols;
 ndim=bf_img0->ndim;
 img0=bf_img0->data;
 img1=bf_img1->data;
 /* Now compute  (r+g+b)/3   r-g  (r+g)/2-b  */

 /* red */
 for(x=0;x<cols;x++) {
   for(y=0;y<rows;y++) {
     cval=img0[y+x*rows];
     img1[y+x*rows]=cval/3.0;
     img1[y+x*rows+rows*cols]=cval;
     img1[y+x*rows+2*rows*cols]=cval/2.0;
   }
 }
 /* green */
 for(x=0;x<cols;x++) {
   for(y=0;y<rows;y++) {
     cval=img0[y+x*rows+rows*cols];
     img1[y+x*rows] +=cval/3.0;
     img1[y+x*rows+rows*cols]+=-cval;
     img1[y+x*rows+2*rows*cols]+=cval/2.0;
   }
 }
 /* blue */
 for(x=0;x<cols;x++) {
   for(y=0;y<rows;y++) {
     cval=img0[y+x*rows+2*rows*cols];
     img1[y+x*rows] +=cval/3.0;
     img1[y+x*rows+2*rows*cols]+=-cval;
   }
 }
}

/*
** Edge list generation function using colour gradient filters
** bfz should be even and >0.
*/
edgeval image_to_edgelist_grad(buffer *bf_image,ebuffer *bf_elist,int bfz) {
int rows,cols,ndim,x,y,m,n,cind,edgeno,nofedges;
int dxz,lpz,dxl,dxh,lpl,lph,cband;
buffer *bf_dx0,*bf_dx1,*bf_dx2;
buffer *bf_derf,*bf_lpf;
fpnum *image,*dx0,*dx1,*dx2,*derf,*lpf;
edge *elist;
fpnum ws;
edgeval d_max;
/* int (*compfn)(const void *,const void *) = NULL;  */
 rows = bf_image->rows; 
 cols = bf_image->cols;
 ndim = bf_image->ndim;
 image = bf_image->data;
 elist = (edge *)bf_elist->data;
 nofedges = bf_elist->cols;
 bf_dx0 = buffer_new(rows,cols,1);
 bf_dx1 = buffer_new(rows,cols,1);
 bf_dx2 = buffer_new(rows,cols,1);
 dx0 = bf_dx0->data;
 dx1 = bf_dx1->data;
 dx2 = bf_dx2->data;

 /* Create filters */
 lpz=bfz-1;dxz=bfz;
 bf_lpf = binfilt1d(lpz);
 lpf = bf_lpf->data;
 bf_derf = binfilt_der(dxz);
 derf = bf_derf->data;

 filter_print(bf_lpf,"lpf");
 filter_print(bf_derf,"derf");
 /* Compute loop ranges for filters */
 dxl=1-dxz/2;dxh=dxl+dxz-1;
 lpl=-(lpz-1)/2;lph=lpl+lpz-1;
 printf("dx_range=[%d,%d]\n",dxl,dxh);
 printf("lp_range=[%d,%d]\n",lpl,lph);
 
 printf("Colour gradients, bfz=%d (should be even)\n",bfz);
  
  d_max=0;  /* Find max */
  for(m=0;m<ndim;m++) {
    cband=m*rows*cols;
    /* Dx for colourband m */
    for(x=0;x<cols-1;x++) {
      for(y=0;y<rows;y++) {
	cind=y+x*rows;
	dx0[cind]=0.0;  /* Overwrite old stuff */
	for(n=dxl;n<=dxh;n++) {
	  if((x+n>=0)&&(x+n<cols)) {
	    dx0[cind]+=derf[n-dxl]*image[cind+rows*n+cband];
	  }
	}
      }
    }
    /* Blur for colourband m */
    for(x=0;x<cols-1;x++) {
      for(y=0;y<rows;y++) {
	cind=y+x*rows;
	dx1[cind]=0.0;   /* Overwrite old stuff */
	ws=0.0;
	for(n=lpl;n<=lph;n++) {
	  if((y+n>=0)&&(y+n<rows)) {
	    dx1[cind]+=lpf[n-lpl]*dx0[cind+n];
	    ws +=lpf[n-lpl];   /* Filter sum */
	  }
	}
	dx1[cind]/=ws;
      }
    }
    /* Square and add to dx2 */
    for(x=0;x<cols-1;x++) {
      for(y=0;y<rows;y++) {
	cind=y+x*rows;
	dx2[cind] += dx1[cind]*dx1[cind];
      }
    }
  }

  /* Build edgelist from Dx */
  edgeno=0;
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      elist[edgeno].dist  = sqrt(dx2[y+x*rows])*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 0;    /* 0 for Dx */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;
      edgeno++;
    }
  }

  /* Clear dx2 */
  for(x=0;x<cols;x++) {
    for(y=0;y<rows;y++) {
      dx2[y+x*rows] = 0.0;
    }
  }

  for(m=0;m<ndim;m++) {
    cband=m*rows*cols;
    /* Dy for colourband m */
    for(x=0;x<cols;x++) {
      for(y=0;y<rows-1;y++) {
	cind=y+x*rows;
	dx0[cind]=0.0;  /* Overwrite old stuff */
	for(n=dxl;n<=dxh;n++) {
	  if((y+n>=0)&&(y+n<rows)) {
	    dx0[cind]+=derf[n-dxl]*image[cind+n+cband];
	  }
	}
      }
    }
    /* Blur for colourband m */
    for(x=0;x<cols;x++) {
      for(y=0;y<rows-1;y++) {
	cind=y+x*rows;
	dx1[cind]=0.0;   /* Overwrite old stuff */
	ws=0.0;
	for(n=lpl;n<=lph;n++) {
	  if((x+n>=0)&&(x+n<cols)) {
	    dx1[cind]+=lpf[n-lpl]*dx0[cind+n*rows];
	    ws +=lpf[n-lpl];   /* Filter sum */
	  }
	}
	dx1[cind]/=ws;
      }
    }
    /* Square and add to dx2 */
    for(x=0;x<cols;x++) {
      for(y=0;y<rows-1;y++) {
	cind=y+x*rows;
	dx2[cind] += dx1[cind]*dx1[cind];
      }
    }
  }

  /* Build edgelist from Dy */
  edgeno=rows*(cols-1);
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = sqrt(dx2[y+x*rows])*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 1;    /* 1 for Dy */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;
      edgeno++;
    }
  }
  /* Free filter buffers */
  buffer_free(bf_dx0);
  buffer_free(bf_dx1);
  buffer_free(bf_dx2);
  /* Free filters */
  buffer_free(bf_derf);
  buffer_free(bf_lpf);
  return d_max;
}
/*
** Edge list generation function using blurred distance maps
*/
edgeval image_to_edgelist_blur(buffer *bf_image,ebuffer *bf_elist,int bfz,int verbosefl) {
int rows,cols,ndim,x,y,cind,edgeno,nofedges;
buffer *bf_dx0;
fpnum *image,*dx0;
edge *elist;
fpnum diff;
edgeval d_max;
/* int (*compfn)(const void *,const void *) = NULL;  */
 rows = bf_image->rows; 
 cols = bf_image->cols;
 ndim = bf_image->ndim;
 image = bf_image->data;
 elist = (edge *)bf_elist->data;
 nofedges = bf_elist->cols;
 bf_dx0 = buffer_new(rows,cols,1);
 dx0 = bf_dx0->data;

 if(verbosefl) printf("Blur distance map, bfz=%d (should be odd)\n",bfz);
  
 d_max=0;  /* Find max */
  /* Loop for Dx */
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      cind=y+x*rows;
      diff=image[cind]-image[cind+rows];                         /* R */
      dx0[cind]  = diff*diff;
      if(ndim>1) {
	diff=image[cind+rows*cols]-image[cind+rows+rows*cols];     /* G */
	dx0[cind] += diff*diff;
	diff=image[cind+2*rows*cols]-image[cind+rows+2*rows*cols]; /* B */
	dx0[cind] += diff*diff;
      }
      dx0[cind] = sqrt(dx0[cind]);                              /* SQRT */
    }
  }

  /* Blur Dx */
  blur_buffer(bf_dx0,bfz);

  /* Build edgelist from Dx */
  edgeno=0;
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 0;    /* 0 for Dx */
     if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }

  /* Loop for Dy (overwrites Dx buffer) */
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      cind=y+x*rows;
      diff=image[cind]-image[cind+1];                         /* R */
      dx0[cind]  = diff*diff;
      if(ndim>1) {
	diff=image[cind+rows*cols]-image[cind+1+rows*cols];     /* G */
	dx0[cind] += diff*diff;
	diff=image[cind+2*rows*cols]-image[cind+1+2*rows*cols]; /* B */
	dx0[cind] += diff*diff;
      }
      dx0[cind] = sqrt(dx0[cind]);                              /* SQRT */
    }
  }
  /* Blur Dy */
  blur_buffer(bf_dx0,bfz);

  /* Build edgelist from Dy */
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 1;    /* 1 for Dy */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }
  /* Free filter buffers */
  buffer_free(bf_dx0);
  return d_max;
}
/*
** Edge list generation function using blurred distance maps
** This version uses 8 neighbours instead of 4.
** Diagonal neighbours are downweighted by 1/sqrt(2)
*/
edgeval image_to_edgelist_blur_n8(buffer *bf_image,ebuffer *bf_elist,int bfz,int verbosefl) {
int rows,cols,ndim,x,y,cind,edgeno,nofedges;
buffer *bf_dx0;
fpnum *image,*dx0;
edge *elist;
fpnum diff;
edgeval d_max;
 rows = bf_image->rows; 
 cols = bf_image->cols;
 ndim = bf_image->ndim;
 image = bf_image->data;
 elist = (edge *)bf_elist->data;
 nofedges = bf_elist->cols;
 bf_dx0 = buffer_new(rows,cols,1);
 dx0 = bf_dx0->data;

 if(nofedges<rows*(cols-1)+cols*(rows-1)+(rows-1)*(cols-1)) {
   printf("Error wrong edgelist length\n");
   return 0;
 }
   
 if(verbosefl) printf("Blur(8) distance map, bfz=%d (should be odd)\n",bfz);
  
 d_max=0;  /* Find max */
  /* Loop for Dx */
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      cind=y+x*rows;
      diff=image[cind]-image[cind+rows];                         /* R */
      dx0[cind]  = diff*diff;
      if(ndim>1) {
	diff=image[cind+rows*cols]-image[cind+rows+rows*cols];     /* G */
	dx0[cind] += diff*diff;
	diff=image[cind+2*rows*cols]-image[cind+rows+2*rows*cols]; /* B */
	dx0[cind] += diff*diff;
      }
      dx0[cind] = sqrt(dx0[cind]);                                 /* SQRT */
    }
  }

  /* Blur Dx */
  blur_buffer(bf_dx0,bfz);

  /* Build edgelist from Dx */
  edgeno=0;
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 0;    /* 0 for Dx */
     if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }

  /* Loop for Dy (overwrite Dx filter buffer) */
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      cind=y+x*rows;
      diff=image[cind]-image[cind+1];                         /* R */
      dx0[cind]  = diff*diff;
      if(ndim>1) {
	diff=image[cind+rows*cols]-image[cind+1+rows*cols];     /* G */
	dx0[cind] += diff*diff;
	diff=image[cind+2*rows*cols]-image[cind+1+2*rows*cols]; /* B */
	dx0[cind] += diff*diff;
      }
      dx0[cind] = sqrt(dx0[cind]);                               /* SQRT */
    }
  }

  /* Blur Dy */
  blur_buffer(bf_dx0,bfz);

  /* Build edgelist from Dy */
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 1;    /* 1 for Dy */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }

  /* Loop for Dxy (overwrite Dy filter buffer) */
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows-1;y++) {
      cind=y+x*rows;
      diff=image[cind]-image[cind+1+rows];                         /* R */
      dx0[cind]  = diff*diff;
      if(ndim>1) {
	diff=image[cind+rows*cols]-image[cind+1+rows+rows*cols];     /* G */
	dx0[cind] += diff*diff;
	diff=image[cind+2*rows*cols]-image[cind+1+rows+2*rows*cols]; /* B */
	dx0[cind] += diff*diff;
      }
      dx0[cind] = sqrt(dx0[cind]/2.0);                         /* SQRT */
    }
  }

  /* Blur Dxy */
  blur_buffer(bf_dx0,bfz);
  /* Build edgelist from Dxy */
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 2;    /* 2 for Dxy */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }

  /* Loop for Dxy2 (overwrite Dxy filter buffer) */
  for(x=1;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      cind=y+x*rows;
      diff=image[cind]-image[cind+1-rows];                         /* R */
      dx0[cind]  = diff*diff;
      if(ndim>1) {
	diff=image[cind+rows*cols]-image[cind+1-rows+rows*cols];     /* G */
	dx0[cind] += diff*diff;
	diff=image[cind+2*rows*cols]-image[cind+1-rows+2*rows*cols]; /* B */
	dx0[cind] += diff*diff;
      }
      dx0[cind] = sqrt(dx0[cind]/2.0);                            /* SQRT */
    }
  }

  /* Blur Dxy */
  blur_buffer(bf_dx0,bfz);
  /* Build edgelist from Dxy */
  for(x=1;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 3;    /* 3 for Dxy2 */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }
/*   printf("edgeno=%d nofedges=%d\n",edgeno,nofedges);  */

  /* Free filter buffer */
  buffer_free(bf_dx0);
  return d_max;
}
/*
** Edge list generation function using blurred distance maps
** This version uses 8 neighbours instead of 4.
** Diagonal neighbours are downweighted by 1/sqrt(2)
** Edges are weighted according to expected variance
*/
edgeval image_to_edgelist_blur_n8_norm(buffer *bf_image,ebuffer *bf_elist,int bfz,int verbosefl) {
int rows,cols,ndim,x,y,cind,edgeno,nofedges;
buffer *bf_dx0;
fpnum *image,*dx0;
edge *elist;
fpnum diff,norm;
edgeval d_max;
/* fpnum wr=1;/\*1.0196e-4;*\/ */
/* fpnum wg=1;/\*4.7059e-5;*\/ */
/* fpnum wb=1;/\*1.5686e-4;*\/ */
/* fpnum q=0;/\*2.5631e-6;   Offset from two quantisation error terms *\/ */

 rows = bf_image->rows; 
 cols = bf_image->cols;
 ndim = bf_image->ndim;
 image = bf_image->data;
 elist = (edge *)bf_elist->data;
 nofedges = bf_elist->cols;
 bf_dx0 = buffer_new(rows,cols,1);
 dx0 = bf_dx0->data;

 if(nofedges<rows*(cols-1)+cols*(rows-1)+(rows-1)*(cols-1)) {
   printf("Error wrong edgelist length\n");
   return 0;
 }
   
 if(verbosefl) printf("Blur(8) normalised distance map, bfz=%d (should be odd)\n",bfz);
  
 d_max=0;  /* Find max */
  /* Loop for Dx */
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      cind=y+x*rows;
      norm=image[cind]+image[cind+rows];                         /* R */
      if(norm>0) {
	diff=image[cind]-image[cind+rows];
	dx0[cind]  = diff*diff/norm;
      } else dx0[cind] = 0;
      if(ndim>1) {
	norm=image[cind+rows*cols]+image[cind+rows+rows*cols];     /* G */
	if(norm>0) {
	  diff=image[cind+rows*cols]-image[cind+rows+rows*cols];
	  dx0[cind] += diff*diff/norm;
	}
	norm=image[cind+2*rows*cols]+image[cind+rows+2*rows*cols]; /* B */
	if(norm>0) {
	  diff=image[cind+2*rows*cols]-image[cind+rows+2*rows*cols];
	  dx0[cind] += diff*diff/norm;
	}
      }
      dx0[cind] = sqrt(dx0[cind]);                               /* SQRT */
    }
  }

  /* Blur Dx */
  blur_buffer(bf_dx0,bfz);

  /* Build edgelist from Dx */
  edgeno=0;
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 0;    /* 0 for Dx */
     if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }

  /* Loop for Dy (overwrite Dx filter buffer) */
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      cind=y+x*rows;
      norm=image[cind]+image[cind+1];                         /* R */
      if(norm>0) {
	diff=image[cind]-image[cind+1];
	dx0[cind] = diff*diff/norm;
      } else dx0[cind] = 0;
      if(ndim>1) {
	norm=image[cind+rows*cols]+image[cind+1+rows*cols];     /* G */
	if(norm>0) {
	  diff=image[cind+rows*cols]-image[cind+1+rows*cols];
	  dx0[cind] += diff*diff/norm;
	}
	norm=image[cind+2*rows*cols]+image[cind+1+2*rows*cols]; /* B */
	if(norm>0) {
	  diff=image[cind+2*rows*cols]-image[cind+1+2*rows*cols];
	  dx0[cind] += diff*diff/norm;
	}
      }
      dx0[cind] = sqrt(dx0[cind]);                            /* SQRT */
    }
  }

  /* Blur Dy */
  blur_buffer(bf_dx0,bfz);

  /* Build edgelist from Dy */
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 1;    /* 1 for Dy */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }

  /* Loop for Dxy (overwrite Dy filter buffer) */
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows-1;y++) {
      cind=y+x*rows;
      norm=image[cind]+image[cind+1+rows];                         /* R */
      if(norm>0) {
	diff=image[cind]-image[cind+1+rows];
	dx0[cind] = diff*diff/norm;
      } else dx0[cind] = 0;
      if(ndim>1) {
	norm=image[cind+rows*cols]+image[cind+1+rows+rows*cols];     /* G */
	if(norm>0) {
	  diff=image[cind+rows*cols]-image[cind+1+rows+rows*cols];
	  dx0[cind] += diff*diff/norm;
	}
	norm=image[cind+2*rows*cols]+image[cind+1+rows+2*rows*cols]; /* B */
	if(norm>0) {
	  diff=image[cind+2*rows*cols]-image[cind+1+rows+2*rows*cols];
	  dx0[cind] += diff*diff/norm;
	}
      }
      dx0[cind] = sqrt(dx0[cind]/2.0);                             /* SQRT */
    }
  }

  /* Blur Dxy */
  blur_buffer(bf_dx0,bfz);
  /* Build edgelist from Dxy */
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 2;    /* 2 for Dxy */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }

  /* Loop for Dxy2 (overwrite Dxy filter buffer) */
  for(x=1;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      cind=y+x*rows;
      norm=image[cind]+image[cind+1-rows];                          /* R */
      if(norm>0) {
	diff=image[cind]-image[cind+1-rows];
	dx0[cind] = diff*diff/norm;
      } else dx0[cind] = 0;
      if(ndim>1) {
	norm=image[cind+rows*cols]+image[cind+1-rows+rows*cols];     /* G */
	if(norm>0) {
	  diff=image[cind+rows*cols]-image[cind+1-rows+rows*cols];
	  dx0[cind] += diff*diff/norm;
	}
	norm=image[cind+2*rows*cols]+image[cind+1-rows+2*rows*cols]; /* B */
	if(norm>0) {
	  diff=image[cind+2*rows*cols]-image[cind+1-rows+2*rows*cols];
	  dx0[cind] += diff*diff/norm;
	}
      }
      dx0[cind] = sqrt(dx0[cind]/2.0);                             /* SQRT */
    }
  }

  /* Blur Dxy */
  blur_buffer(bf_dx0,bfz);
  /* Build edgelist from Dxy */
  for(x=1;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 3;    /* 3 for Dxy2 */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }
/*   printf("edgeno=%d nofedges=%d\n",edgeno,nofedges);  */

  /* Free filter buffer */
  buffer_free(bf_dx0);
  return d_max;
}

/*
** Edge list generation function using blurred distance maps
** Edges are weighted according to expected variance
*/
edgeval image_to_edgelist_blur_norm(buffer *bf_image,ebuffer *bf_elist,int bfz,int verbosefl) {
int rows,cols,ndim,x,y,cind,edgeno,nofedges;
buffer *bf_dx0,*bf_dx1;
fpnum *image,*dx0,*dx1;
edge *elist;
fpnum diff,norm;
edgeval d_max;
/* fpnum wr=1;/\*1.0196e-4;*\/ */
/* fpnum wg=1;/\*4.7059e-5;*\/ */
/* fpnum wb=1;/\*1.5686e-4;*\/ */
/* fpnum q=0;/\*2.5631e-6;   Offset from two quantisation error terms *\/ */

 rows = bf_image->rows; 
 cols = bf_image->cols;
 ndim = bf_image->ndim;
 image = bf_image->data;
 elist = (edge *)bf_elist->data;
 nofedges = bf_elist->cols;
 bf_dx0 = buffer_new(rows,cols,1);
 bf_dx1 = buffer_new(rows,cols,1);
 dx0 = bf_dx0->data;
 dx1 = bf_dx1->data;

 if(verbosefl) printf("Blur normalised distance map, bfz=%d (should be odd)\n",bfz);
  
 d_max=0;  /* Find max */
  /* Loop for Dx */
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      cind=y+x*rows;
      norm=image[cind]+image[cind+rows];                         /* R */
      if(norm>0) {
	diff=image[cind]-image[cind+rows];
	dx0[cind] += diff*diff/norm;
      } else dx0[cind] = 0;
      if(ndim>1) {
	norm=image[cind+rows*cols]+image[cind+rows+rows*cols];     /* G */
	if(norm>0) {
	  diff=image[cind+rows*cols]-image[cind+rows+rows*cols];
	  dx0[cind] += diff*diff/norm;
	}
	norm=image[cind+2*rows*cols]+image[cind+rows+2*rows*cols]; /* B */
	if(norm>0) {
	  diff=image[cind+2*rows*cols]-image[cind+rows+2*rows*cols];
	  dx0[cind] += diff*diff/norm;
	}
      }
      dx0[cind] = sqrt(dx0[cind]);                               /* SQRT */
    }
  }

  /* Blur Dx */
  blur_buffer(bf_dx0,bfz);

  /* Build edgelist from Dx */
  edgeno=0;
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 0;    /* 0 for Dx */
     if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }

  /* Loop for Dy */
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      cind=y+x*rows;
      norm=image[cind]+image[cind+1];                         /* R */
      if(norm>0) {
	diff=image[cind]-image[cind+1];
	dx1[cind]  = diff*diff/norm;
      } else dx1[cind] = 0;
      if(ndim>1) {
	norm=image[cind+rows*cols]+image[cind+1+rows*cols];     /* G */
	if(norm>0) {
	  diff=image[cind+rows*cols]-image[cind+1+rows*cols];
	  dx1[cind] += diff*diff/norm;
	}
	norm=image[cind+2*rows*cols]+image[cind+1+2*rows*cols]; /* B */
	if(norm>0) {
	  diff=image[cind+2*rows*cols]-image[cind+1+2*rows*cols];
	  dx1[cind] += diff*diff/norm;
	}
      }
      dx1[cind] = sqrt(dx1[cind]);                            /* SQRT */
    }
  }
  /* Blur Dy */
  blur_buffer(bf_dx1,bfz);

  /* Build edgelist from Dy */
  edgeno=rows*(cols-1);
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = dx1[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 1;    /* 1 for Dy */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }
  /* Free filter buffers */
  buffer_free(bf_dx0);
  buffer_free(bf_dx1);
  return d_max;
}
/*
** Edge list generation function using blurred distance maps
*/
edgeval image_to_edgelist_blur8(buffer *bf_image,buffer *bf_elist,int bfz) {
int rows,cols,ndim,x,y,cind,edgeno,nofedges;
buffer *bf_dx0,*bf_dx1;
unsigned char *image;
fpnum *dx0,*dx1;
edge *elist;
fpnum diff;
edgeval d_max;
/* int (*compfn)(const void *,const void *) = NULL;  */
 rows = bf_image->rows; 
 cols = bf_image->cols;
 ndim = bf_image->ndim;
 image = (unsigned char *)bf_image->data;
 elist = (edge *)bf_elist->data;
 nofedges = bf_elist->cols;
 bf_dx0 = buffer_new(rows,cols,1);
 bf_dx1 = buffer_new(rows,cols,1);
 dx0 = bf_dx0->data;
 dx1 = bf_dx1->data;

 printf("image_to_edgelist_blur: bfz=%d (should be odd)\n",bfz);
  
 d_max=0;  /* Find max */
  /* Loop for Dx */
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      cind=y+x*rows;
      diff=image[cind]-image[cind+rows];                         /* R */
      dx0[cind]  = diff*diff;
      diff=image[cind+rows*cols]-image[cind+rows+rows*cols];     /* G */
      dx0[cind] += diff*diff;
      diff=image[cind+2*rows*cols]-image[cind+rows+2*rows*cols]; /* B */
      dx0[cind] = sqrt(dx0[cind]+diff*diff)/255.0;               /* SQRT */
    }
  }

  /* Blur Dx */
  blur_buffer(bf_dx0,bfz);

  /* Build edgelist from Dx */
  edgeno=0;
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      elist[edgeno].dist  = dx0[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 0;    /* 0 for Dx */
     if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }

  /* Loop for Dy */
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      cind=y+x*rows;
      diff=image[cind]-image[cind+1];                         /* R */
      dx1[cind]  = diff*diff;
      diff=image[cind+rows*cols]-image[cind+1+rows*cols];     /* G */
      dx1[cind] += diff*diff;
      diff=image[cind+2*rows*cols]-image[cind+1+2*rows*cols]; /* B */
      dx1[cind] = sqrt(dx1[cind]+diff*diff)/255.0;            /* SQRT */
    }
  }
  /* Blur Dy */
  blur_buffer(bf_dx1,bfz);

  /* Build edgelist from Dy */
  edgeno=rows*(cols-1);
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = dx1[y+x*rows]*EDGE_SCALE+EDGE_OFFSET;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 1;    /* 1 for Dy */
      if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist;     
      edgeno++;
    }
  }
  /* Free filter buffers */
  buffer_free(bf_dx0);
  buffer_free(bf_dx1);
  return d_max;
}
/*
** Old edge list generation function (Squared distances)
*/
fpnum image_to_edgelist0(buffer *bf_image,buffer *bf_elist) {
int rows,cols,ndim,x,y,m,edgeno,nofedges;
fpnum *image;
fpnum diff,d_max;
edge *elist;
 rows = bf_image->rows; 
 cols = bf_image->cols;
 ndim = bf_image->ndim;
 image = bf_image->data;
 elist = (edge *)bf_elist->data;
 nofedges = bf_elist->cols;

  printf("image_to_edgelist.\n");

  /* First loop for Dx */
  edgeno=0;
  for(x=0;x<cols-1;x++) {
    for(y=0;y<rows;y++) {
      diff=image[y+x*rows]-image[y+(x+1)*rows];
      elist[edgeno].dist  = diff*diff;
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 0;    /* 0 for Dx */
      edgeno++;
    }
  }
  /* Remaining loops for Dx */
  for(m=1;m<ndim;m++) {
    edgeno=0;
    for(x=0;x<cols-1;x++) {
      for(y=0;y<rows;y++) {
	diff = image[y+x*rows+m*rows*cols]-image[y+(x+1)*rows+m*rows*cols];
	elist[edgeno].dist += diff*diff;
	edgeno++;
      }
    }
  }
  /* First loop for Dy */
  edgeno=rows*(cols-1);
  for(x=0;x<cols;x++) {
    for(y=0;y<rows-1;y++) {
      elist[edgeno].dist  = SQR(image[y+x*rows]-image[y+1+x*rows]);
      elist[edgeno].xpos  = x+1;  /* Matlab indices start at 1 */
      elist[edgeno].ypos  = y+1;
      elist[edgeno].dirfl = 1;    /* 1 for Dy */
      edgeno++;
    }
  }
  /* Remaining loops for Dy */
  for(m=1;m<ndim;m++) {
    edgeno=rows*(cols-1);
    for(x=0;x<cols;x++) {
      for(y=0;y<rows-1;y++) {
	diff = image[y+x*rows+m*rows*cols]-image[y+1+x*rows+m*rows*cols];
	elist[edgeno].dist += diff*diff;
	edgeno++;
      }
    }
  }
  /* Round, and find max */
  d_max=0;
  for(edgeno=0;edgeno<nofedges;edgeno++) {
    elist[edgeno].dist=floor(elist[edgeno].dist*1e10+0.5)/1e10;
    if(elist[edgeno].dist>d_max) d_max=elist[edgeno].dist; /* NB! Squared */
  }
  return d_max;
}
/*
** Edge list update function
*/
void update_edgelist(buffer *bf_image,edge *elist,int nedge,
		     ibuffer *bf_limg,buffer *bf_pvec,buffer *bf_mvec) {
int *limg;
int edgeno,lno,ind1,ind2,x,y,rows,cols;
fpnum *pvec,*image,*mvec;
fpnum pvec1[3],pvec2[3];
fpnum dist,diff,a;
/*   printf("update_edgelist\n"); */
  limg = bf_limg->data;
  rows = bf_limg->rows;
  cols = bf_limg->cols;
  pvec = bf_pvec->data;
  mvec = bf_mvec->data;
  image = bf_image->data;
  for(edgeno=0;edgeno<nedge;edgeno++) {
    x=elist[edgeno].xpos;
    y=elist[edgeno].ypos;
    if(elist[edgeno].dirfl==0) {
      /* Dx */
      ind1=y-1+(x-1)*rows;
      ind2=y-1+x*rows;
    } else {
      /* Dy */
      ind1=y-1+(x-1)*rows; 
      ind2=y+(x-1)*rows; 
    }
    /* Pos 1 */
    lno=limg[ind1];
    if(lno==0) {
      /* No region, read from image */
      pvec1[0]=image[ind1];
      pvec1[1]=image[ind1+rows*cols];
      pvec1[2]=image[ind1+rows*cols*2];
    } else {
      a=mvec[lno*6+0];
      pvec1[0]=pvec[lno*3+0]/a;
      pvec1[1]=pvec[lno*3+1]/a;
      pvec1[2]=pvec[lno*3+2]/a;
    }
    /* Pos 2 */
    lno=limg[ind2];
    if(lno==0) {
      /* No region, read from image */
      pvec2[0]=image[ind2];
      pvec2[1]=image[ind2+rows*cols];
      pvec2[2]=image[ind2+rows*cols*2];
    } else {
      a=mvec[lno*6+0];
      pvec2[0]=pvec[lno*3+0]/a;
      pvec2[1]=pvec[lno*3+1]/a;
      pvec2[2]=pvec[lno*3+2]/a;
    }
    /* Compute new distance */
    diff=pvec1[0]-pvec2[0];dist  = diff*diff;
    diff=pvec1[1]-pvec2[1];dist += diff*diff;
    diff=pvec1[2]-pvec2[2];dist += diff*diff;
    elist[edgeno].dist = dist;    /* Squared distance */
  }
}
/*
 *  Remove holes in mvec,bbox and dalist
 *  and update labels in limg
 */
void bloblist_defrag(buffer *bf_mvec,buffer *bf_pvec,ibuffer *bf_bbox,
		     buffer *bf_dalist,ibuffer *bf_limg,int *pnbl,int *plno) {
  ibuffer *bf_indl;
  int k,rind;
  int x,y,rows,cols,cind;
  int *limg,*indl;
  fpnum *pvec;
  da_tuple *dalist;
  boundingbox *bbox;
  momentvector *mvec;
  mvec=(momentvector *)bf_mvec->data;
  pvec=bf_pvec->data;
  limg=bf_limg->data;
  bbox=(boundingbox *)bf_bbox->data;
  dalist=(da_tuple *)bf_dalist->data;
  rows=bf_limg->rows;
  cols=bf_limg->cols;
  bf_indl=ibuffer_new(plno[0],1,1);
  indl=bf_indl->data;
  rind=1;  /* Start at 1, label 0 means unassigned */
  /* Generate index list and compact mvec and bbox */
  indl[0]=0;
  for(k=1;k<plno[0];k++) {
    if(mvec[k].area>0) {
      if(k>rind) {
	mvec[rind]   = mvec[k];
	bbox[rind]   = bbox[k];
	dalist[rind] = dalist[k];
	memcpy(&pvec[rind*3],&pvec[k*3],3*sizeof(fpnum));
      }
      indl[k]=rind++;
    }
  }
  /* Relabel limg */
  for(x=0;x<cols;x++) {
    for(y=0;y<rows;y++) {
      cind=y+rows*x;
      limg[cind]=indl[limg[cind]];
    }
  }
  pnbl[0]=0;     /* No holes left */
  plno[0]=rind;  /* New next label */
  ibuffer_free(bf_indl);
}
/*
** Region growing evolution function (defrag)
*/
void edgelist_to_bloblist00(buffer **bfp_mvec,buffer **bfp_pvec,buffer **bfp_arate,buffer *bf_image,buffer *bf_elist,buffer *bf_thres,int amin,double size_inc,double defrag_frac,int defrag_min,fpnum res)
{
int rows,cols,first_edge,edges_left,cnedge;
int edgeno,k,boutcnt,ind1,ind2,cx,cy,xx,yy,ccx,ccy;
int max_label,lno,l1,l2,cl,lh,lno_old;
fpnum *image;
edge *elist; 
fpnum d_pivot,new_arate,old_arate;
int *limg;
buffer *bf_mvec,*bf_aold,*bf_mvec_out,*bf_dalist,*bf_arate; 
buffer *bf_pvec,*bf_pvec_out;
fpnum *aold,*thres=NULL;
momentvector *mvec,*mvec_out; 
arearate *arate;
da_tuple *dalist;
fpnum *pvec_out,*pvec,*pvec_k,*pvec_l;
ibuffer *bf_limg,*bf_bbox;
boundingbox *bbox; 
int timestep,nbl;
int xmin,xmax,ymin,ymax,out_ind;
int tslist_flag;

  printf("Defrag version...\n");
  tslist_flag=(bf_thres!=NULL);  /* Did we get a threshold list? */
  if(tslist_flag) {
    thres=bf_thres->data;
    printf("Using adaptive timestep\n");
  } else {
    printf("Using fixed timestep res=%g\n",res);
  }
  rows=bf_image->rows;
  cols=bf_image->cols; 
  image=bf_image->data;
  bf_limg=ibuffer_new(rows,cols,1);    /* Label image */
  limg=bf_limg->data;
  elist=(edge *)bf_elist->data;
  edges_left=bf_elist->cols;
  first_edge=0;
  printf("edges_left=%d\n",edges_left);
  printf("amin=%d size_inc=%g\n",amin,size_inc);
  printf("defrag_frac=%g defrag_min=%d\n",defrag_frac,defrag_min);
  /* Look at input image size */
  buffer_pdims(bf_image);
  /* Number of table entries */
  max_label=(rows*cols*3)/10;
  printf("max_label=%d\n",max_label);
  /* Mvec and Pvec arrays */
  bf_mvec_out=buffer_new(6,max_label,1);
  mvec_out=(momentvector *)bf_mvec_out->data;
  bf_pvec_out=buffer_new(3,max_label,1);
  pvec_out=bf_pvec_out->data;
  bf_mvec=buffer_new(6,max_label,1);
  mvec=(momentvector *)bf_mvec->data;
  bf_pvec=buffer_new(3,max_label,1);
  pvec=bf_pvec->data;
  /* Other arrays */
  bf_aold=buffer_new(1,max_label,1);   /* Old areas */
  aold=bf_aold->data;
  bf_dalist=buffer_new(3,max_label,1);     /* rows: dist0,area0,int:(ts0,oind) */
  dalist=(da_tuple *)bf_dalist->data;
  bf_arate=buffer_new(6,max_label,1);      /* rows: (arate,amin,d0,dn,t0,tn) */
  arate=(arearate *)bf_arate->data;
  bf_bbox=ibuffer_new(4,max_label,1);      /* rows: xmin,xmax,ymin,ymax */
  bbox=(boundingbox *)bf_bbox->data;

  lno=1;        /* Running label number */
  lno_old=lno;  /* Highest label in last evolution step */
  boutcnt=0;    /* Number of blobs stored */
  nbl=0;        /* Count number of blanks (for defrag) */

  d_pivot=0.0;  /* Current distance */
  timestep=0;   /* Count evolution steps */

  while(edges_left) {
    /* Move small valued edges to front */
    cnedge=0;
    while(cnedge==0) {
      if(tslist_flag) {
/* 	d_pivot=sqrt(thres[timestep]); */
	d_pivot=thres[timestep];
      } else {
	d_pivot+=res;
      }
      cnedge=edgelist_partition((edge *)&elist[first_edge],
				edges_left,d_pivot);
/* 				edges_left,SQR(d_pivot)); */
      timestep++;
    }
    /* Go through all edges in this time step */
    for(edgeno=first_edge;edgeno<first_edge+cnedge;edgeno++) {
      cx=elist[edgeno].xpos-1;  /* Remove matlab offset */
      cy=elist[edgeno].ypos-1;
      ind1=cy+rows*cx;
      if(elist[edgeno].dirfl==0) {
	/* Dx */
	ind2=ind1+rows;
      } else {
	/* Dy */
	ind2=ind1+1;
      }
      /* Read labels */
      l1=limg[ind1];
      l2=limg[ind2];
      if((l1==0)||(l2==0)) {
	if((l1==0)&&(l2==0)) {
	  /* New region */
	  cl=lno++;
	  if(lno>max_label) {lno=max_label; printf("OOPS!\n");}
	  /* Remember at what threshold we appeared */
	  dalist[cl].dist = d_pivot;
	  dalist[cl].area = 2;
	  dalist[cl].ts   = timestep;
	  dalist[cl].oind = 0;     /* Region not stored yet */
	  /* Set bounding box, mvec and pvec */
	  pvec_k = &pvec[cl*3];
	  ccx=cx+1;ccy=cy+1;   /* Add offset for moment computation */
	  if(elist[edgeno].dirfl==0) {
	    /* Dx */
	    bbox[cl].xmin=cx;
	    bbox[cl].xmax=cx+1;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+(ccx+1);
	    mvec[cl].my=ccy+ccy;
	    mvec[cl].mx2=2*ccx*(ccx+1)+1;  /* cx*cx+(cx+1)*(cx+1) */
	    mvec[cl].mxy=2*ccx*ccy+ccy;    /* cx*cy+(cx+1)*cy */
	    mvec[cl].my2=2*ccy*ccy;
	  } else {
	    /* Dy */
	    bbox[cl].xmin=cx;
	    bbox[cl].xmax=cx;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy+1;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+ccx;
	    mvec[cl].my=ccy+(ccy+1);
	    mvec[cl].mx2=2*ccx*ccx; 
	    mvec[cl].mxy=2*ccx*ccy+ccx;     /* cx*cy+cx*(cy+1) */
	    mvec[cl].my2=2*ccy*(ccy+1)+1;   /* cy*cy+(cy+1)*(cy+1) */
	  }
	  /* Initialize pvec */
	  pvec_k[0]=image[ind1]+image[ind2];
	  pvec_k[1]=image[ind1+rows*cols]+image[ind2+rows*cols];
	  pvec_k[2]=image[ind1+rows*cols*2]+image[ind2+rows*cols*2];
	  /* Set labels */
	  limg[ind1]=cl;
	  limg[ind2]=cl;      
	} else {
	  /* Append one pixel to a region */
	  if(l1>l2) {cl=l1;} else {cl=l2;} /* cl=max(l1,l2), non-zero label */
	  ccx=cx+1;ccy=cy+1;   /* Add offset for moment computation */
	  /* Update bounding box */
	  if(elist[edgeno].dirfl==0) {
	    /* Dx */
	    if(bbox[cl].xmin>cx)   bbox[cl].xmin=cx;
	    if(bbox[cl].xmax<cx+1) bbox[cl].xmax=cx+1;
	    if(cl==l1) ccx++;   /* Right one is new (cx+1,cy) */
	  } else {
	    /* Dy */
	    if(bbox[cl].ymin>cy)   bbox[cl].ymin=cy;
	    if(bbox[cl].ymax<cy+1) bbox[cl].ymax=cy+1;
	    if(cl==l1) ccy++;   /* Lower one is new (cx,cy+1) */
	  }
	  /* Update mvec */
	  mvec[cl].area += 1;
	  mvec[cl].mx   += ccx;
	  mvec[cl].my   += ccy;
	  mvec[cl].mx2  += ccx*ccx;
	  mvec[cl].mxy  += ccx*ccy;
	  mvec[cl].my2  += ccy*ccy;
	  
	  ind1=ccy-1+(ccx-1)*rows;   /* Index of new pixel */
	  /* Update pvec */
	  pvec_k = &pvec[cl*3];
	  pvec_k[0] += image[ind1];
	  pvec_k[1] += image[ind1+rows*cols];
	  pvec_k[2] += image[ind1+rows*cols*2];
	  /* Set label */
	  limg[ind1] = cl;
	}
      } else {
	if(l1 != l2) {
	  /* Merge regions */
	  /* Largest region keeps label */
 	  if(mvec[l1].area>mvec[l2].area) {cl=l1;lh=l2;} else {cl=l2;lh=l1;} 
	  /* Change all lh to cl */
	  pvec_k = &pvec[lh*3];
	  xmin=bbox[lh].xmin;
	  xmax=bbox[lh].xmax;
	  ymin=bbox[lh].ymin;
	  ymax=bbox[lh].ymax;
	  for(xx=xmin;xx<=xmax;xx++) {
	    for(yy=ymin;yy<=ymax;yy++) {
	      if(limg[yy+xx*rows]==lh) limg[yy+xx*rows]=cl;
	    }
	  }
	  /* Merge bounding boxes */
	  if(bbox[cl].xmin>bbox[lh].xmin) bbox[cl].xmin=bbox[lh].xmin;
	  if(bbox[cl].xmax<bbox[lh].xmax) bbox[cl].xmax=bbox[lh].xmax;
	  if(bbox[cl].ymin>bbox[lh].ymin) bbox[cl].ymin=bbox[lh].ymin;
	  if(bbox[cl].ymax<bbox[lh].ymax) bbox[cl].ymax=bbox[lh].ymax;
	  /* Merge mvecs */
	  mvec[cl].area += mvec[lh].area;
	  mvec[cl].mx   += mvec[lh].mx;
	  mvec[cl].my   += mvec[lh].my;
	  mvec[cl].mx2  += mvec[lh].mx2;
	  mvec[cl].mxy  += mvec[lh].mxy;
	  mvec[cl].my2  += mvec[lh].my2;
	  /* Merge pvecs */
	  pvec_l = &pvec[cl*3];
	  pvec_l[0] += pvec_k[0];
	  pvec_l[1] += pvec_k[1];
	  pvec_l[2] += pvec_k[2];
	  /* Indicate as merged */
	  mvec[lh].area=0;
	  nbl++;  /* One more blank */
	}
      }
    }
    /*
     * New step in evolution, check all regions
     */
/*     cl=limg[285+215*rows]; */
/*     new_arate=(mvec[cl].area-dalist[cl].area)/(d_pivot-dalist[cl].dist); */
/*     printf("%d: dalist=[%g %g %d] nfl=%d\tarate=%g \ta=%g ao=%g cfl=%d\n",timestep,dalist[cl].area,dalist[cl].dist,dalist[cl].ts,mvec[cl].area>aold[cl]*size_inc,new_arate,mvec[cl].area,aold[cl],mvec[cl].area<dalist[cl].area*2.0); */
    for(k=0;k<lno_old;k++) {
      /* Require min area */
      if(mvec[k].area>amin) {
	/* Require small area increase since start */
	if(mvec[k].area<dalist[k].area*2.0) {
/*   	  new_arate=mvec[k].area/dalist[k].area/(d_pivot-dalist[k].dist);  */
	  if(timestep-dalist[k].ts>2) {
	    new_arate=(mvec[k].area-dalist[k].area)/(d_pivot-dalist[k].dist); 
	  } else {
	    new_arate=1e14;
	  }
	  /* Find out where to store region */
	  out_ind=dalist[k].oind;
	  if(!out_ind) {
	    /* New out_ind */
	    out_ind=boutcnt++;
	    dalist[k].oind=out_ind;
	    /* Remember start posn */
	    arate[out_ind].d0 = dalist[k].dist; /* d_0 */
	    arate[out_ind].t0 = dalist[k].ts;   /* t_0 */
	    old_arate=1e15;
	  } else {
	    old_arate=arate[out_ind].arate;
	  }
	  if(new_arate<old_arate) { 
	    /* Update blob description */
	    mvec_out[out_ind]=mvec[k];
	    memcpy(&pvec_out[out_ind*3],&pvec[k*3],3*sizeof(fpnum));
	    /* Update opt point */
	    arate[out_ind].arate = new_arate;
	    arate[out_ind].amin  = mvec[k].area;
 	  }
	  /* Update end point */
	  arate[out_ind].dn    = d_pivot;
 	  arate[out_ind].tn    = timestep; 
	} else {
	  if(dalist[k].oind) {
	    if(mvec[k].area<aold[k]*size_inc) {
	      mvec_out[dalist[k].oind].area=0;         /* Remove region */
	    }
	  }
	}
	if(boutcnt>max_label) {
	  printf("OOPS!");
	  boutcnt=max_label;
	}
      }
      /* Check if large size increase */
      if(mvec[k].area>aold[k]*size_inc) {
	/* NOTE: We cannot set end point here, this will miss mergers */
	dalist[k].dist = d_pivot;      /* Reset start distance */
	dalist[k].area = mvec[k].area; /* Reset start area */
	dalist[k].ts   = timestep;     /* Reset start timestep */
	dalist[k].oind = 0;            /* Terminate */
      }
    }
    /* Run defrag if necessary */
    if((nbl>lno*defrag_frac)&&(nbl>defrag_min)) {
      bloblist_defrag(bf_mvec,bf_pvec,bf_bbox,bf_dalist,bf_limg,&nbl,&lno);
    }
    /* Copy mvec to aold */
    for(k=0;k<lno;k++) {
      aold[k]=mvec[k].area;
    }
    
    /* End of while loop */
    lno_old=lno;
    first_edge += cnedge;
    edges_left -= cnedge;
      
/*      printf("%d ",cnedge);  */
/*     if(cnedge>0) { */
/*       update_edgelist(bf_image,(edge *)&(elist[first_edge].dist), */
/* 		      edges_left,bf_limg,bf_pvec,bf_mvec); */
/*     } */
  }
  ibuffer_free(bf_bbox);
  ibuffer_free(bf_limg);
  printf("boutcnt=%d\n",boutcnt);
  printf("timestep=%d d_pivot=%g\n",timestep,d_pivot);
  buffer_free(bf_mvec);
  buffer_free(bf_aold);
  buffer_free(bf_dalist);
  buffer_free(bf_pvec);
  bfp_mvec[0]=bf_mvec_out;
  bfp_pvec[0]=bf_pvec_out;
  bfp_arate[0]=bf_arate;
}
/*
 * Path compression recursion
 */
int find_next(int *next,int k) {
  int retval;
  retval=k;
  if(next[k] != k) {
    retval=find_next(next,next[k]);
    next[k]=retval;
  }
  return retval;
}
/*
** Region growing evolution function (forward ptrs and path compression)
*/
void edgelist_to_bloblist0(buffer **bfp_mvec,buffer **bfp_pvec,buffer **bfp_arate,buffer *bf_image,buffer *bf_elist,buffer *bf_thres,int amin,double size_inc,double defrag_frac,int defrag_min,fpnum res)
{
int rows,cols,first_edge,edges_left,cnedge;
int edgeno,k,boutcnt,ind1,ind2,cx,cy,xx,yy,ccx,ccy;
int max_label,l1,l2,cl,lh,free_label;
fpnum *image;
edge *elist; 
fpnum d_pivot,new_arate,old_arate;
int *limg;
buffer *bf_mvec,*bf_aold,*bf_mvec_out,*bf_dalist,*bf_arate; 
buffer *bf_pvec,*bf_pvec_out;
fpnum *aold,*thres=NULL;
momentvector *mvec,*mvec_out; 
arearate *arate;
da_tuple *dalist;
fpnum *pvec_out,*pvec,*pvec_k,*pvec_l;
ibuffer *bf_limg,*bf_bbox;
ibuffer *bf_next;
int *next;
boundingbox *bbox; 
int timestep;
int xmin,xmax,ymin,ymax,out_ind;
int tslist_flag;

  printf("Path compression version...\n");
  tslist_flag=(bf_thres!=NULL);  /* Did we get a threshold list? */
  if(tslist_flag) {
    thres=bf_thres->data;
    printf("Using adaptive timestep\n");
  } else {
    printf("Using fixed timestep res=%g\n",res);
  }
  rows=bf_image->rows;
  cols=bf_image->cols; 
  image=bf_image->data;
  bf_limg=ibuffer_new(rows,cols,1);    /* Label image */
  limg=bf_limg->data;
  elist=(edge *)bf_elist->data;
  edges_left=bf_elist->cols;
  first_edge=0;
  printf("edges_left=%d\n",edges_left);
  printf("amin=%d size_inc=%g\n",amin,size_inc);
  printf("defrag_frac=%g defrag_min=%d\n",defrag_frac,defrag_min);
  /* Look at input image size */
  buffer_pdims(bf_image);
  /* Number of table entries */
  max_label=(rows*cols*3)/10;
  /* Mvec and Pvec arrays */
  bf_mvec_out=buffer_new(6,max_label,1);
  mvec_out=(momentvector *)bf_mvec_out->data;
  bf_pvec_out=buffer_new(3,max_label,1);
  pvec_out=bf_pvec_out->data;
  bf_mvec=buffer_new(6,max_label,1);
  mvec=(momentvector *)bf_mvec->data;
  bf_pvec=buffer_new(3,max_label,1);
  pvec=bf_pvec->data;
  /* Other arrays */
  bf_aold=buffer_new(1,max_label,1);   /* Old areas */
  aold=bf_aold->data;
  bf_dalist=buffer_new(3,max_label,1);     /* rows: dist0,area0,int:(ts0,oind) */
  dalist=(da_tuple *)bf_dalist->data;
  bf_arate=buffer_new(6,max_label,1);      /* rows: (arate,amin,d0,dn,t0,tn) */
  arate=(arearate *)bf_arate->data;
  bf_bbox=ibuffer_new(4,max_label,1);      /* rows: xmin,xmax,ymin,ymax */
  bbox=(boundingbox *)bf_bbox->data;
  bf_next=ibuffer_new(1,max_label,1);      /* rows: next */
  next=(int *)bf_next->data;

  for(k=0;k<max_label;k++) next[k]=k;      /* Initialise next pointers */

  free_label=1; /* Running label number */
  boutcnt=0;    /* Number of blobs stored */

  d_pivot=0.0;  /* Current distance */
  timestep=0;   /* Count evolution steps */

  while(edges_left) {
    /* Move small valued edges to front */
    cnedge=0;
    while(cnedge==0) {
      if(tslist_flag) {
/* 	d_pivot=sqrt(thres[timestep]); */
	d_pivot=thres[timestep];
      } else {
	d_pivot+=res;
      }
      cnedge=edgelist_partition((edge *)&elist[first_edge],
				edges_left,d_pivot);
/* 				edges_left,SQR(d_pivot)); */
      timestep++;
    }
    /* Go through all edges in this time step */
    for(edgeno=first_edge;edgeno<first_edge+cnedge;edgeno++) {
      cx=elist[edgeno].xpos-1;  /* Remove matlab offset */
      cy=elist[edgeno].ypos-1;
      ind1=cy+rows*cx;
      if(elist[edgeno].dirfl==0) {
	/* Dx */
	ind2=ind1+rows;
      } else {
	/* Dy */
	ind2=ind1+1;
      }
      /* Read labels */
      l1=limg[ind1];
      l2=limg[ind2];
      if((l1==0)||(l2==0)) {
	if((l1==0)&&(l2==0)) {
	  /* New region */
	  cl=free_label++;    /* New label */

	  if(free_label>max_label) {free_label=max_label; printf("OOPS!\n");}
	  /* Remember at what threshold we appeared */
	  dalist[cl].dist = d_pivot;
	  dalist[cl].area = 2;
	  dalist[cl].ts   = timestep;
	  dalist[cl].oind = 0;     /* Region not stored yet */
	  /* Set bounding box, mvec and pvec */
	  pvec_k = &pvec[cl*3];
	  ccx=cx+1;ccy=cy+1;   /* Add offset for moment computation */
	  if(elist[edgeno].dirfl==0) {
	    /* Dx */
	    bbox[cl].xmin=cx;
	    bbox[cl].xmax=cx+1;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+(ccx+1);
	    mvec[cl].my=ccy+ccy;
	    mvec[cl].mx2=2*ccx*(ccx+1)+1;  /* cx*cx+(cx+1)*(cx+1) */
	    mvec[cl].mxy=2*ccx*ccy+ccy;    /* cx*cy+(cx+1)*cy */
	    mvec[cl].my2=2*ccy*ccy;
	  } else {
	    /* Dy */
	    bbox[cl].xmin=cx;
	    bbox[cl].xmax=cx;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy+1;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+ccx;
	    mvec[cl].my=ccy+(ccy+1);
	    mvec[cl].mx2=2*ccx*ccx; 
	    mvec[cl].mxy=2*ccx*ccy+ccx;     /* cx*cy+cx*(cy+1) */
	    mvec[cl].my2=2*ccy*(ccy+1)+1;   /* cy*cy+(cy+1)*(cy+1) */
	  }
	  /* Initialize pvec */
	  pvec_k[0]=image[ind1]+image[ind2];
	  pvec_k[1]=image[ind1+rows*cols]+image[ind2+rows*cols];
	  pvec_k[2]=image[ind1+rows*cols*2]+image[ind2+rows*cols*2];
	  /* Set labels */
	  limg[ind1]=cl;
	  limg[ind2]=cl;      
	} else {
	  /* Append one pixel to a region */
	  if(l1>l2) {cl=l1;} else {cl=l2;} /* cl=max(l1,l2), non-zero label */
	  ccx=cx+1;ccy=cy+1;   /* Add offset for moment computation */
	  /* Update bounding box */
	  if(elist[edgeno].dirfl==0) {
	    /* Dx */
	    if(bbox[cl].xmin>cx)   bbox[cl].xmin=cx;
	    if(bbox[cl].xmax<cx+1) bbox[cl].xmax=cx+1;
	    if(cl==l1) ccx++;   /* Right one is new (cx+1,cy) */
	  } else {
	    /* Dy */
	    if(bbox[cl].ymin>cy)   bbox[cl].ymin=cy;
	    if(bbox[cl].ymax<cy+1) bbox[cl].ymax=cy+1;
	    if(cl==l1) ccy++;   /* Lower one is new (cx,cy+1) */
	  }
	  /* Update mvec */
	  mvec[cl].area += 1;
	  mvec[cl].mx   += ccx;
	  mvec[cl].my   += ccy;
	  mvec[cl].mx2  += ccx*ccx;
	  mvec[cl].mxy  += ccx*ccy;
	  mvec[cl].my2  += ccy*ccy;
	  
	  ind1=ccy-1+(ccx-1)*rows;   /* Index of new pixel */
	  /* Update pvec */
	  pvec_k = &pvec[cl*3];
	  pvec_k[0] += image[ind1];
	  pvec_k[1] += image[ind1+rows*cols];
	  pvec_k[2] += image[ind1+rows*cols*2];
	  /* Set label */
	  limg[ind1] = cl;
	}
      } else {
	if(l1 != l2) {
	  /* Merge regions */
	  /* Largest region keeps label */
 	  if(mvec[l1].area>mvec[l2].area) {cl=l1;lh=l2;} else {cl=l2;lh=l1;} 
	  /* Change all lh to cl */
	  pvec_k = &pvec[lh*3];
	  xmin=bbox[lh].xmin;
	  xmax=bbox[lh].xmax;
	  ymin=bbox[lh].ymin;
	  ymax=bbox[lh].ymax;
	  for(xx=xmin;xx<=xmax;xx++) {
	    for(yy=ymin;yy<=ymax;yy++) {
	      if(limg[yy+xx*rows]==lh) limg[yy+xx*rows]=cl;
	    }
	  }
	  /* Merge bounding boxes */
	  if(bbox[cl].xmin>bbox[lh].xmin) bbox[cl].xmin=bbox[lh].xmin;
	  if(bbox[cl].xmax<bbox[lh].xmax) bbox[cl].xmax=bbox[lh].xmax;
	  if(bbox[cl].ymin>bbox[lh].ymin) bbox[cl].ymin=bbox[lh].ymin;
	  if(bbox[cl].ymax<bbox[lh].ymax) bbox[cl].ymax=bbox[lh].ymax;
	  /* Merge mvecs */
	  mvec[cl].area += mvec[lh].area;
	  mvec[cl].mx   += mvec[lh].mx;
	  mvec[cl].my   += mvec[lh].my;
	  mvec[cl].mx2  += mvec[lh].mx2;
	  mvec[cl].mxy  += mvec[lh].mxy;
	  mvec[cl].my2  += mvec[lh].my2;
	  /* Merge pvecs */
	  pvec_l = &pvec[cl*3];
	  pvec_l[0] += pvec_k[0];
	  pvec_l[1] += pvec_k[1];
	  pvec_l[2] += pvec_k[2];
	  /* Indicate as merged */
	  mvec[lh].area=0;
	  next[lh]=next[lh+1];  /* Point to what next is pointing to */
	}
      }
    }
    /*
     * New step in evolution, check all regions
     */
/*     cl=limg[285+215*rows]; */
/*     new_arate=(mvec[cl].area-dalist[cl].area)/(d_pivot-dalist[cl].dist); */
/*     printf("%d: dalist=[%g %g %d] nfl=%d\tarate=%g \ta=%g ao=%g cfl=%d\n",timestep,dalist[cl].area,dalist[cl].dist,dalist[cl].ts,mvec[cl].area>aold[cl]*size_inc,new_arate,mvec[cl].area,aold[cl],mvec[cl].area<dalist[cl].area*2.0); */
    for(k=find_next(next,1);k<free_label;k=find_next(next,k+1)) {

      /* Require min area */
      if(mvec[k].area>amin) {
	/* Require small area increase since start */
	if(mvec[k].area<dalist[k].area*2.0) {
/*   	  new_arate=mvec[k].area/dalist[k].area/(d_pivot-dalist[k].dist);  */
	  if(timestep-dalist[k].ts>2) {
	    new_arate=(mvec[k].area-dalist[k].area)/(d_pivot-dalist[k].dist); 
	  } else {
	    new_arate=1e14;
	  }
	  /* Find out where to store region */
	  out_ind=dalist[k].oind;
	  if(!out_ind) {
	    /* New out_ind */
	    out_ind=boutcnt++;
	    dalist[k].oind=out_ind;
	    /* Remember start posn */
	    arate[out_ind].d0 = dalist[k].dist; /* d_0 */
	    arate[out_ind].t0 = dalist[k].ts;   /* t_0 */
	    old_arate=1e15;
	  } else {
	    old_arate=arate[out_ind].arate;
	  }
	  if(new_arate<old_arate) { 
	    /* Update blob description */
	    mvec_out[out_ind]=mvec[k];
	    memcpy(&pvec_out[out_ind*3],&pvec[k*3],3*sizeof(fpnum));
	    /* Update opt point */
	    arate[out_ind].arate = new_arate;
	    arate[out_ind].amin  = mvec[k].area;
 	  }
	  /* Update end point */
	  arate[out_ind].dn    = d_pivot;
 	  arate[out_ind].tn    = timestep; 
	} else {
	  if(dalist[k].oind) {
	    if(mvec[k].area<aold[k]*size_inc) {
	      mvec_out[dalist[k].oind].area=0;         /* Remove region */
	    }
	  }
	}
	if(boutcnt>max_label) {
	  printf("OOPS!");
	  boutcnt=max_label;
	}
      }
      /* Check if large size increase */
      if(mvec[k].area>aold[k]*size_inc) {
	/* NOTE: We cannot set end point here, this will miss mergers */
	/* NOTE: Regions start with area=2, then get correct area here */
	dalist[k].dist = d_pivot;      /* Reset start distance */
	dalist[k].area = mvec[k].area; /* Reset start area */
	dalist[k].ts   = timestep;     /* Reset start timestep */
	dalist[k].oind = 0;            /* Terminate */
      }
      aold[k]=mvec[k].area; /* Set aold */
    }
    /* End of while loop */
    first_edge += cnedge;
    edges_left -= cnedge;      
  }

  ibuffer_free(bf_bbox);
  ibuffer_free(bf_limg);
  ibuffer_free(bf_next);
  printf("boutcnt=%d\n",boutcnt);
  printf("free_label=%d max_label=%d\n",free_label,max_label);
  printf("timestep=%d d_pivot=%g\n",timestep,d_pivot);
  buffer_free(bf_mvec);
  buffer_free(bf_aold);
  buffer_free(bf_dalist);
  buffer_free(bf_pvec);
  bfp_mvec[0]=bf_mvec_out;
  bfp_pvec[0]=bf_pvec_out;
  bfp_arate[0]=bf_arate;
}
/*
** Region growing evolution function (doubly linked list)
*/
void edgelist_to_bloblist(buffer **bfp_mvec,buffer **bfp_pvec,buffer **bfp_arate,buffer *bf_image,ebuffer *bf_elist,ebuffer *bf_thres,int amin,double size_inc,fpnum res,int verbosefl)
{
int rows,cols,ndim,first_edge,edges_left,cnedge;
int edgeno,k,boutcnt,ind1,ind2,cx,cy,xx,yy,ccx,ccy;
int max_label,l1,l2,cl,lh,free_label;
fpnum *image;
/* unsigned char *image; */
edge *elist; 
fpnum new_arate,old_arate;
int *limg;
buffer *bf_mvec,*bf_aold,*bf_mvec_out,*bf_dalist,*bf_arate; 
buffer *bf_pvec,*bf_pvec_out;
fpnum *aold;
edgeval *thres=NULL,d_pivot;
momentvector *mvec,*mvec_out; 
arearate *arate;
da_tuple *dalist;
fpnum *pvec_out,*pvec,*pvec_k,*pvec_l;
ibuffer *bf_limg,*bf_bbox;
ibuffer *bf_nodes;
node *active_list,*nodes,*cnode;
boundingbox *bbox; 
int timestep;
int xmin,xmax,ymin,ymax,out_ind;
int tslist_flag;

  if(verbosefl) printf("Doubly-linked-list version...\n");
  tslist_flag=(bf_thres!=NULL);  /* Did we get a threshold list? */
  if(tslist_flag) {
    thres=bf_thres->data;
    if(verbosefl) printf("Using adaptive timestep\n");
  } else {
    if(verbosefl) printf("Using fixed timestep res=%g\n",res);
  }
  rows=bf_image->rows;
  cols=bf_image->cols; 
  ndim=bf_image->ndim; 
  image=bf_image->data;
  bf_limg=ibuffer_new(rows,cols,1);    /* Label image */
  limg=bf_limg->data;
  elist=(edge *)bf_elist->data;
  edges_left=bf_elist->cols;
  first_edge=0;
  /* Number of table entries */
  max_label=rows*cols;/*(rows*cols*3)/10;*/
  /* Mvec and Pvec arrays */
  bf_mvec_out=buffer_new(6,max_label,1);
  mvec_out=(momentvector *)bf_mvec_out->data;
  bf_pvec_out=buffer_new(3,max_label,1);
  pvec_out=bf_pvec_out->data;
  bf_mvec=buffer_new(6,max_label,1);
  mvec=(momentvector *)bf_mvec->data;
  bf_pvec=buffer_new(3,max_label,1);
  pvec=bf_pvec->data;
  /* Other arrays */
  bf_aold=buffer_new(1,max_label,1);   /* Old areas */
  aold=bf_aold->data;
  bf_dalist=buffer_new(3,max_label,1);     /* rows: dist0,area0,int:(ts0,oind) */
  dalist=(da_tuple *)bf_dalist->data;
  bf_arate=buffer_new(6,max_label,1);      /* rows: (arate,amin,d0,dn,t0,tn) */
  arate=(arearate *)bf_arate->data;
  bf_bbox=ibuffer_new(4,max_label,1);      /* rows: xmin,xmax,ymin,ymax */
  bbox=(boundingbox *)bf_bbox->data;
  bf_nodes=ibuffer_new(2,max_label,1);     /* rows: prev,next */
  nodes=(node *)bf_nodes->data;
  active_list=NULL;
/*   active_end=NULL; */

  free_label=1; /* Running label number */
  boutcnt=0;    /* Number of blobs stored */

  d_pivot=0.0;  /* Current distance */
  timestep=0;   /* Count evolution steps */

  while(edges_left) {
    /* Move small valued edges to front */
    cnedge=0;
    while(cnedge==0) {
      if(tslist_flag) {
/* 	d_pivot=sqrt(thres[timestep]); */
	d_pivot=thres[timestep];
      } else {
	d_pivot+=res;
      }
      cnedge=edgelist_partition(&elist[first_edge],edges_left,d_pivot);
/* 				edges_left,SQR(d_pivot)); */
      timestep++;
    }
    /* Go through all edges in this time step */
    for(edgeno=first_edge;edgeno<first_edge+cnedge;edgeno++) {
      cx=elist[edgeno].xpos-1;  /* Remove matlab offset */
      cy=elist[edgeno].ypos-1;
      ind1=cy+rows*cx;
      switch((int)elist[edgeno].dirfl) {
      case 0:	/* Dx */
	ind2=ind1+rows;
	break;
      case 1:   /* Dy */	
	ind2=ind1+1;
	break;
      case 2:   /* Dxy */	
	ind2=ind1+1+rows;
	break;
      default:/* Dxy2 */	
	ind2=ind1+1-rows;
	break;
      }
      /* Read labels */
      l1=limg[ind1];
      l2=limg[ind2];
      if(l1==l2) {
	if(l1==0) {
	  /* New region */
	  cl=free_label++;    /* New label */

	  /* Insert at end of list */
/* 	  nodes[cl].next=NULL; */
/* 	  nodes[cl].prev=active_end; */
/* 	  if(active_end != NULL) active_end->next = &nodes[cl]; */
/* 	  active_end = &nodes[cl]; */
/* 	  if(active_list == NULL) active_list= &nodes[cl]; */

	  /* Insert at start of list */
	  nodes[cl].next = active_list;
	  nodes[cl].prev = NULL;
	  if(active_list != NULL) active_list->prev = &nodes[cl];
	  active_list = &nodes[cl];

	  if(free_label>max_label) {free_label=max_label; printf("OOPS!\n");}
	  /* Remember at what threshold we appeared */
	  dalist[cl].dist = d_pivot;
	  dalist[cl].area = 2;
	  dalist[cl].ts   = timestep;
	  dalist[cl].oind = 0;     /* Region not stored yet */
	  /* Set bounding box, mvec and pvec */
	  pvec_k = &pvec[cl*3];
	  ccx=cx+1;ccy=cy+1;   /* Add offset for moment computation */
	  switch((int)elist[edgeno].dirfl) {
	  case 0:  /* Dx */
	    bbox[cl].xmin=cx;
	    bbox[cl].xmax=cx+1;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+(ccx+1);
	    mvec[cl].my=ccy+ccy;
	    mvec[cl].mx2=2*ccx*(ccx+1)+1;  /* cx*cx+(cx+1)*(cx+1) */
	    mvec[cl].mxy=2*ccx*ccy+ccy;    /* cx*cy+(cx+1)*cy */
	    mvec[cl].my2=2*ccy*ccy;
	    break;
	  case 1:  /* Dy */
	    bbox[cl].xmin=cx;
	    bbox[cl].xmax=cx;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy+1;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+ccx;
	    mvec[cl].my=ccy+(ccy+1);
	    mvec[cl].mx2=2*ccx*ccx; 
	    mvec[cl].mxy=2*ccx*ccy+ccx;     /* cx*cy+cx*(cy+1) */
	    mvec[cl].my2=2*ccy*(ccy+1)+1;   /* cy*cy+(cy+1)*(cy+1) */
	    break;
	  case 2:  /* Dxy */
	    bbox[cl].xmin=cx;
	    bbox[cl].xmax=cx+1;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy+1;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+(ccx+1);
	    mvec[cl].my=ccy+(ccy+1);
	    mvec[cl].mx2=2*ccx*(ccx+1)+1;    /* cx*cx+(cx+1)*(cx+1) */
	    mvec[cl].mxy=2*ccx*ccy+ccx+ccy+1;/* cx*cy+(cx+1)*(cy+1) */
	    mvec[cl].my2=2*ccy*(ccy+1)+1;    /* cy*cy+(cy+1)*(cy+1) */
	    break;
	  default: /* Dxy2 */
	    bbox[cl].xmin=cx-1;
	    bbox[cl].xmax=cx;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy+1;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+(ccx-1);
	    mvec[cl].my=ccy+(ccy+1);
	    mvec[cl].mx2=2*ccx*(ccx-1)+1;    /* (cx-1)*(cx-1)+cx*cx */
	    mvec[cl].mxy=2*ccx*ccy+ccx-ccy-1;/* cx*cy+(cx-1)*(cy+1) */
	    mvec[cl].my2=2*ccy*(ccy+1)+1;    /* cy*cy+(cy+1)*(cy+1) */
	  }
	  /* Initialize pvec */
	  pvec_k[0]=image[ind1]+image[ind2];
	  if(ndim>1) {
	    pvec_k[1]=image[ind1+rows*cols]+image[ind2+rows*cols];
	    pvec_k[2]=image[ind1+rows*cols*2]+image[ind2+rows*cols*2];
	  }
	  /* Set labels */
	  limg[ind1]=cl;
	  limg[ind2]=cl;      
	}
      } else {
	if((l1==0)||(l2==0)) {
	  /* Append one pixel to a region */
	  if(l1>l2) {cl=l1;} else {cl=l2;} /* cl=max(l1,l2), non-zero label */
	  ccx=cx+1;ccy=cy+1;   /* Add offset for moment computation */
	  /* Update bounding box */
	  switch((int)elist[edgeno].dirfl) {
	  case 0:    /* Dx */
	    if(bbox[cl].xmin>cx)   bbox[cl].xmin=cx;
	    if(bbox[cl].xmax<cx+1) bbox[cl].xmax=cx+1;
	    if(cl==l1) ccx++;   /* Right one is new (cx+1,cy) */
	    break;
	  case 1:    /* Dy */
	    if(bbox[cl].ymin>cy)   bbox[cl].ymin=cy;
	    if(bbox[cl].ymax<cy+1) bbox[cl].ymax=cy+1;
	    if(cl==l1) ccy++;   /* Lower one is new (cx,cy+1) */
	    break;
	  case 2:    /* Dxy */
	    if(bbox[cl].ymin>cy)   bbox[cl].ymin=cy;
	    if(bbox[cl].ymax<cy+1) bbox[cl].ymax=cy+1;
	    if(bbox[cl].xmin>cx)   bbox[cl].xmin=cx;
	    if(bbox[cl].xmax<cx+1) bbox[cl].xmax=cx+1;
	    if(cl==l1) {ccx++;ccy++;}  /* Lower one is new (cx+1,cy+1) */
	    break;
	  default:   /* Dxy2 */
	    if(bbox[cl].ymin>cy)   bbox[cl].ymin=cy;
	    if(bbox[cl].ymax<cy+1) bbox[cl].ymax=cy+1;
	    if(bbox[cl].xmin>cx-1) bbox[cl].xmin=cx-1;
	    if(bbox[cl].xmax<cx)   bbox[cl].xmax=cx;
	    if(cl==l1) {ccx--;ccy++;}  /* Lower one is new (cx-1,cy+1) */
	  }
	  /* Update mvec */
	  mvec[cl].area += 1;
	  mvec[cl].mx   += ccx;
	  mvec[cl].my   += ccy;
	  mvec[cl].mx2  += ccx*ccx;
	  mvec[cl].mxy  += ccx*ccy;
	  mvec[cl].my2  += ccy*ccy;
	  
	  ind1=ccy-1+(ccx-1)*rows;   /* Index of new pixel */
	  /* Update pvec */
	  pvec_k = &pvec[cl*3];
	  pvec_k[0] += image[ind1];
	  if(ndim>1) {
	    pvec_k[1] += image[ind1+rows*cols];
	    pvec_k[2] += image[ind1+rows*cols*2];
	  }
	  /* Set label */
	  limg[ind1] = cl;
	} else {
	  /* Merge regions */
	  /* Largest region keeps label */
 	  if(mvec[l1].area>mvec[l2].area) {cl=l1;lh=l2;} else {cl=l2;lh=l1;} 
	  /* Change all lh to cl */
	  pvec_k = &pvec[lh*3];
	  xmin=bbox[lh].xmin;
	  xmax=bbox[lh].xmax;
	  ymin=bbox[lh].ymin;
	  ymax=bbox[lh].ymax;
	  for(xx=xmin;xx<=xmax;xx++) {
	    for(yy=ymin;yy<=ymax;yy++) {
	      if(limg[yy+xx*rows]==lh) limg[yy+xx*rows]=cl;
	    }
	  }
	  /* Merge bounding boxes */
	  if(bbox[cl].xmin>bbox[lh].xmin) bbox[cl].xmin=bbox[lh].xmin;
	  if(bbox[cl].xmax<bbox[lh].xmax) bbox[cl].xmax=bbox[lh].xmax;
	  if(bbox[cl].ymin>bbox[lh].ymin) bbox[cl].ymin=bbox[lh].ymin;
	  if(bbox[cl].ymax<bbox[lh].ymax) bbox[cl].ymax=bbox[lh].ymax;
	  /* Merge mvecs */
	  mvec[cl].area += mvec[lh].area;
	  mvec[cl].mx   += mvec[lh].mx;
	  mvec[cl].my   += mvec[lh].my;
	  mvec[cl].mx2  += mvec[lh].mx2;
	  mvec[cl].mxy  += mvec[lh].mxy;
	  mvec[cl].my2  += mvec[lh].my2;
	  /* Merge pvecs */
	  pvec_l = &pvec[cl*3];
	  pvec_l[0] += pvec_k[0];
	  pvec_l[1] += pvec_k[1];
	  pvec_l[2] += pvec_k[2];
	  /* Indicate as merged */
	  mvec[lh].area=0;

	  /* Unlink node from list */
	  if(nodes[lh].prev != NULL) 
	    nodes[lh].prev->next = nodes[lh].next;
	  else
	    active_list = nodes[lh].next; /* First in list removed */
	  if(nodes[lh].next != NULL) 
	    nodes[lh].next->prev = nodes[lh].prev;
/* 	  else */
/* 	    active_end = nodes[lh].prev; /\* Last in list removed *\/ */
	}
      }
    }
    /*
     * New step in evolution, check all regions
     */
/*     cl=limg[285+215*rows]; */
/*     new_arate=(mvec[cl].area-dalist[cl].area)/(d_pivot-dalist[cl].dist); */
/*     printf("%d: dalist=[%g %g %d] nfl=%d\tarate=%g \ta=%g ao=%g cfl=%d\n",timestep,dalist[cl].area,dalist[cl].dist,dalist[cl].ts,mvec[cl].area>aold[cl]*size_inc,new_arate,mvec[cl].area,aold[cl],mvec[cl].area<dalist[cl].area*2.0); */
/*     for(k=0;k<top_label;k++) { */
    cnode = active_list;
    while(cnode != NULL) {
      /* Get index */
      k = cnode-nodes;
      cnode = cnode->next; /* Next node to visit */
      /* Require min area */
      if(mvec[k].area>amin) {
	/* Require small area increase since start */
	if(mvec[k].area<dalist[k].area*2.0) {
/*   	  new_arate=mvec[k].area/dalist[k].area/(d_pivot-dalist[k].dist);  */
	  if(timestep-dalist[k].ts>2) {
	    new_arate=(mvec[k].area-dalist[k].area)/(d_pivot-dalist[k].dist); 
	  } else {
	    new_arate=1e14;
	  }
	  /* Find out where to store region */
	  out_ind=dalist[k].oind;
	  if(!out_ind) {
	    /* New out_ind */
	    out_ind=boutcnt++;
	    dalist[k].oind=out_ind;
	    /* Remember start posn */
	    arate[out_ind].d0 = dalist[k].dist; /* d_0 */
	    arate[out_ind].t0 = dalist[k].ts;   /* t_0 */
	    old_arate=1e15;
	  } else {
	    old_arate=arate[out_ind].arate;
	  }
	  if(new_arate<old_arate) { 
	    /* Update blob description */
	    mvec_out[out_ind]=mvec[k];
	    memcpy(&pvec_out[out_ind*3],&pvec[k*3],3*sizeof(fpnum));
	    /* Update opt point */
	    arate[out_ind].arate = new_arate;
	    arate[out_ind].amin  = mvec[k].area;
 	  }
	  /* Update end point */
	  arate[out_ind].dn    = d_pivot;
 	  arate[out_ind].tn    = timestep; 
	} else {
	  if(dalist[k].oind) {
	    if(mvec[k].area<aold[k]*size_inc) {
	      mvec_out[dalist[k].oind].area=0;         /* Remove region */
	    }
	  }
	}
	if(boutcnt>max_label) {
	  printf("OOPS!");
	  boutcnt=max_label;
	}
      }
      /* Check if large size increase */
      if(mvec[k].area>aold[k]*size_inc) {
	/* NOTE: We cannot set end point here, this will miss mergers */
	/* NOTE: Regions start with area=2, then get correct area here */
	dalist[k].dist = d_pivot;      /* Reset start distance */
	dalist[k].area = mvec[k].area; /* Reset start area */
	dalist[k].ts   = timestep;     /* Reset start timestep */
	dalist[k].oind = 0;            /* Terminate */
      }
      aold[k]=mvec[k].area; /* Set aold */
    }
    /* End of while loop */
    first_edge += cnedge;
    edges_left -= cnedge;
      
/*      printf("%d ",cnedge);  */
/*     if(cnedge>0) { */
/*       update_edgelist(bf_image,&elist[first_edge], */
/* 		      edges_left,bf_limg,bf_pvec,bf_mvec); */
/*     } */
  }
  ibuffer_free(bf_bbox);
  ibuffer_free(bf_limg);
  ibuffer_free(bf_nodes);
/*   printf("boutcnt=%d\n",boutcnt); */
/*   printf("timestep=%d d_pivot=%g\n",timestep,d_pivot); */
  buffer_free(bf_mvec);
  buffer_free(bf_aold);
  buffer_free(bf_dalist);
  buffer_free(bf_pvec);
  bfp_mvec[0]=bf_mvec_out;
  bfp_pvec[0]=bf_pvec_out;
  bfp_arate[0]=bf_arate;
}

/*
** Distance cdf generation function
*/
void edgelist_to_cdf(buffer *bf_elist,buffer *bf_ecdf,edgeval d_max) {
int Ns,k,l,nedges;
fpnum *ecdf;
edge *elist;
edgeval cdist; 
  Ns     = bf_ecdf->cols;
  ecdf   = bf_ecdf->data;
  nedges = bf_elist->cols;
  elist  = (edge *)bf_elist->data;
  printf("edgelist_to_cdf: d_max=%g Ns=%d\n",(fpnum)d_max,Ns);
  /* Compute list of d sample points */
  for(k=1;k<Ns-1;k++) {
    ecdf[k*2]=d_max*exp(log(2.0)*(k-Ns+1-1));
  }
  ecdf[(Ns-1)*2]=d_max;  /* Set last sample point to .5*2*d_max */

  /* Sample corresponding cdf values */
  for(k=0;k<nedges;k++) {
    cdist=elist[k].dist;
    for(l=1;l<Ns-1;l++) {
      if(cdist<ecdf[l*2]) ecdf[l*2+1]+=1.0;
    }
  }
  /* Normalise */
  for(l=1;l<Ns-1;l++) {
    ecdf[l*2+1]/=nedges;
  }
  ecdf[0+1]=0.0;          /* Set first value */
  ecdf[(Ns-1)*2+1]=1.0;   /* Set last value */
}
/*
** Evolution thresholds computation function
*  (Old, uses samples on log scale)
*/
void evolution_thresholds0(buffer *bf_ecdf,buffer *bf_thres) {
int N,Ns,k,cind;
fpnum cc=0.0;
fpnum *ecdf,*thres;
  N     = bf_thres->cols; /* Nof thresholds to compute */
  Ns    = bf_ecdf->cols;
  ecdf  = bf_ecdf->data;
  thres = bf_thres->data;
  printf("Old evolution (using sampled cdf)\n");
  cind=0;  /* Index of right value in interpolation */
  for(k=0;k<N-1;k++) {
    cc=(fpnum)(k+1.0)/(fpnum)N; /* current fractile */
    while(cc>ecdf[cind*2+1]) cind++;
    thres[k] =exp(((cc-ecdf[(cind-1)*2+1])*log(ecdf[cind*2])+
	       (ecdf[cind*2+1]-cc)*log(ecdf[(cind-1)*2]))/
             (ecdf[cind*2+1]-ecdf[(cind-1)*2+1]));
  }
  thres[N-1]=4*ecdf[2*(Ns-1)];  /* Set last value */
  printf("cind=%d Ns=%d\n",cind,Ns);
  printf("cc=%g ecdf[cind*2+1]=%g\n",cc,ecdf[cind*2+1]);
}
/*
** Evolution thresholds computation function
*  uses Exponential distr.
*/
double evolution_thresholds(ebuffer *bf_elist,ebuffer *bf_thres) {
int N,Ne,k;
fpnum cc;
double emean;
edgeval *thres;
edge *elist;
  N     = bf_thres->cols; /* Nof thresholds to compute */
  Ne    = bf_elist->cols;
  elist = (edge *)bf_elist->data;
  thres = (edgeval *)bf_thres->data;
  emean=0.0;
  for(k=0;k<Ne;k++) emean += elist[k].dist;
  emean /= (double)Ne;
  for(k=0;k<N-1;k++) {
    cc=(fpnum)(k+1.0)/(fpnum)N; /* current fractile */
/*     thres[k] = sqrt(-log(1-cc)*emean*emean*4.0*PI); Rayleigh cdf */ 
    thres[k] =-log(1-cc)*emean;   /* Exponential cdf */
  }
#ifdef INTEGER_EDGES
  thres[N-1]=INT_MAX;
#else
  thres[N-1]=100*thres[N-2];     /* Last thres should be = inf */
#endif
  return emean;
}
/*
** Evolution thresholds computation function
*  uses Chi^2 distr.
*/
double evolution_thresholds2(ebuffer *bf_elist,ebuffer *bf_thres,int order) {
int N,Ne,k,tind;
fpnum cc;
double emean,dist;
edgeval *thres;
edge *elist;
#include "chi_table.h"
  N     = bf_thres->cols; /* Nof thresholds to compute */
  Ne    = bf_elist->cols;
  elist = (edge *)bf_elist->data;
  thres = (edgeval *)bf_thres->data;
  emean=0.0;
  for(k=0;k<Ne;k++) emean += elist[k].dist;
  emean /= (double)Ne;
  switch(order) {
  case 1:
    /* chitab1[0] = cinv(0), chitab1[TABLE_SIZE-1] = cinv(1) */
    for(k=0;k<N-1;k++) {
      cc=(fpnum)(k+1.0)/(fpnum)N;    /* Current fractile */
      tind=floor(cc*(TABLE_SIZE-1)); /* Closest index before */
      dist=cc*(TABLE_SIZE-1)-tind;   /* Remainder */
      thres[k]=emean*(chitab1[tind]*(1-dist)+chitab1[tind+1]*dist);
    }
    break;
  case 2:
    for(k=0;k<N-1;k++) {
      cc=(fpnum)(k+1.0)/(fpnum)N;    /* Current fractile */
      thres[k] =-log(1-cc)*emean;    /* Exponential cdf (Chi^2_2) */
    }
    break;
  case 3:
    /* chitab3[0] = cinv(0), chitab3[TABLE_SIZE-1] = cinv(1) */
    for(k=0;k<N-1;k++) {
      cc=(fpnum)(k+1.0)/(fpnum)N;    /* Current fractile */
      tind=floor(cc*(TABLE_SIZE-1)); /* Closest index before */
      dist=cc*(TABLE_SIZE-1)-tind;   /* Remainder */
      thres[k]=emean*(chitab3[tind]*(1-dist)+chitab3[tind+1]*dist);
    }
  }
#ifdef INTEGER_EDGES
  thres[N-1]=INT_MAX;
#else
  thres[N-1]=100*thres[N-2];     /* Last thres should be = inf */
#endif
  return emean;
}
/*
 *  Destructively convert raw moments to centered moments
 *
 *  bf_mvec     Moment vectors
 *  bf_pvec     Colour (property) moment vectors
 *
 */
void center_moments(buffer *bf_mvec,buffer *bf_pvec) {
  fpnum *pvec,*pvec_k;
  momentvector *mvec;
  fpnum a,mx,my;
  int k,blobs;
  mvec=(momentvector *)bf_mvec->data;
  pvec=bf_pvec->data;
  blobs=bf_mvec->cols;
  /* Convert raw moments to centroid, inertia etc. */
  for(k=0;k<blobs;k++) {
    a=mvec[k].area;
    if(a>0) {
      mvec[k].mx/=a;mx=mvec[k].mx;
      mvec[k].my/=a;my=mvec[k].my;
      mvec[k].mx2=mvec[k].mx2/a-mx*mx;
      mvec[k].mxy=mvec[k].mxy/a-mx*my;
      mvec[k].my2=mvec[k].my2/a-my*my;
      pvec_k=&pvec[k*3];
      pvec_k[0]/=a;
      pvec_k[1]/=a;
      pvec_k[2]/=a;
    }
  }
}
/*
 *  Require that all blobs have:
 *  1. det(I)>0, 2. mvec(0)>amin, 3. margin>min_margin
 *
 *  Mark invalid blobs by setting mvec(0)=0
 *  and count number of valid blobs.
 *
 *  bf_mvec     Moment vectors
 *  amin        Minimum required area
 *  bf_arate    Arate vectors
 *  min_margin  Minimum required margin
 */
int bloblist_mark_invalid(buffer *bf_mvec,int amin,buffer *bf_arate,
			  fpnum min_margin) {
  int regions,k,validcnt;
  arearate *arate;
  momentvector *mvec;
  mvec    = (momentvector *)bf_mvec->data;
  regions = bf_mvec->cols;
  arate   = (arearate *)bf_arate->data;

  validcnt=0;
  for(k=0;k<regions;k++) {
    if((mvec[k].area>amin)&&((mvec[k].mx2*mvec[k].my2-SQR(mvec[k].mxy))>2e-9)
       &&((arate[k].dn-arate[k].d0)>min_margin)) {
      validcnt++;
    } else {
      mvec[k].area=0; /* Mark as invalid */
    }
  }
  return(validcnt);
}
/*
 *  Require that all blobs have a minor axis longer than 1.5 pixels
 *
 *  Mark invalid blobs by setting mvec(0)=0
 *  and count number of valid blobs.
 *
 *  bf_mvec     Moment vectors
 */
int bloblist_shape_invalid(buffer *bf_mvec) {
  int regions,k,validcnt;
  fpnum a,b,lambda2;
  momentvector *mvec;
  mvec    = (momentvector *)bf_mvec->data;
  regions = bf_mvec->cols;

  validcnt=0;
  for(k=0;k<regions;k++) {
    if(mvec[k].area>0) {
      /* Find smallest eigenvalue */ 
      a=mvec[k].mx2+mvec[k].my2;                         /* trace */
      b=sqrt(SQR(mvec[k].mx2-mvec[k].my2)+4*SQR(mvec[k].mxy));  
/*      lambda1=(a+b)/2; */
      lambda2=(a-b)/2;
/*       if(2*sqrt(lambda2)>1.5) */
      if(lambda2>0.5625) validcnt++;
      else mvec[k].area=0; /* Mark as invalid */
    }
  }
  return(validcnt);
}
/*
 * Copy valid blobs (blobs with mvec(0)>0) 
 * from bf_mvec to bf_mvec2
 * and from bf_arate to bf_arate2
 *
 */
void bloblist_compact(buffer *bf_mvec,buffer *bf_mvec2,
		      buffer *bf_pvec,buffer *bf_pvec2,
		      buffer *bf_arate,buffer *bf_arate2) {
  fpnum *pvec,*pvec2;
  momentvector *mvec,*mvec2;
  arearate *arate,*arate2;
  int rind,blobs,k;
  pvec   = bf_pvec->data;
  pvec2  = bf_pvec2->data;
  mvec   = (momentvector *)bf_mvec->data;
  mvec2  = (momentvector *)bf_mvec2->data;
  arate  = (arearate *)bf_arate->data;
  arate2 = (arearate *)bf_arate2->data;
  blobs  = bf_mvec->cols;
  rind=0;
  for(k=0;k<blobs;k++) {
    if(mvec[k].area>0) {
      mvec2[rind]  = mvec[k];
      arate2[rind] = arate[k];
      memcpy(&pvec2[rind*3],&pvec[k*3],3*sizeof(fpnum));
      rind++;
    }
  }
}
/*
 * Copy valid blobs (blobs with mvec(0)>0) 
 * from bf_mvec to bf_mvec2
 * and from bf_pvec to bf_pvec2
 *
 */
void bloblist_compact2(buffer *bf_mvec,buffer *bf_mvec2,
		       buffer *bf_pvec,buffer *bf_pvec2) {
  fpnum *pvec,*pvec2;
  momentvector *mvec,*mvec2;
  int rind,blobs,k;
  pvec   = bf_pvec->data;
  pvec2  = bf_pvec2->data;
  mvec   = (momentvector *)bf_mvec->data;
  mvec2  = (momentvector *)bf_mvec2->data;
  blobs  = bf_mvec->cols;
  rind=0;
  for(k=0;k<blobs;k++) {
    if(mvec[k].area>0) {
      mvec2[rind]  = mvec[k];
      memcpy(&pvec2[rind*3],&pvec[k*3],3*sizeof(fpnum));
      rind++;
    }
  }
}



/*
** Region growing evolution function (doubly linked list) with mask foreground/baskground 
 * Michela Farenzena, September 2009
*/
void edgelist_to_bloblist_masked(buffer **bfp_mvec,buffer **bfp_pvec,buffer **bfp_arate,buffer *bf_image, buffer *bf_mask, ebuffer *bf_elist,ebuffer *bf_thres,int amin,double size_inc,fpnum res,int verbosefl)
{

    int rows,cols,ndim,first_edge,edges_left,cnedge;
int edgeno,k,boutcnt,ind1,ind2,cx,cy,xx,yy,ccx,ccy;
int max_label,l1,l2,cl,lh,free_label;
fpnum *image; fpnum *mask;

/* unsigned char *image; */
edge *elist; 
fpnum new_arate,old_arate;
int *limg;
buffer *bf_mvec,*bf_aold,*bf_mvec_out,*bf_dalist,*bf_arate; 
buffer *bf_pvec,*bf_pvec_out;
fpnum *aold;
edgeval *thres=NULL,d_pivot;
momentvector *mvec,*mvec_out; 
arearate *arate;
da_tuple *dalist;
fpnum *pvec_out,*pvec,*pvec_k,*pvec_l;
ibuffer *bf_limg,*bf_bbox;
ibuffer *bf_nodes;
node *active_list,*nodes,*cnode;
boundingbox *bbox; 
int timestep;
int xmin,xmax,ymin,ymax,out_ind;
int tslist_flag;

  if(verbosefl) printf("Doubly-linked-list version...\n");
  tslist_flag=(bf_thres!=NULL);  /* Did we get a threshold list? */
  if(tslist_flag) {
    thres=bf_thres->data;
    if(verbosefl) printf("Using adaptive timestep\n");
  } else {
    if(verbosefl) printf("Using fixed timestep res=%g\n",res);
  }
  
  rows=bf_image->rows;
  cols=bf_image->cols; 
  ndim=bf_image->ndim; 
  image=bf_image->data;
  mask = bf_mask->data;
  
  bf_limg=ibuffer_new(rows,cols,1);    /* Label image */
  limg=bf_limg->data;
  
  // Set Label(i,j) = -1 if Mask(i,j) == 0 (background pixel
  for(k=0; k< rows*cols; k++)
  { if ( mask[k] == 0 )
        limg[k] = -1 ;   
  }
  
  elist=(edge *)bf_elist->data;
  edges_left=bf_elist->cols;
  first_edge=0;
  
/* Number of table entries */
  max_label=rows*cols;/*(rows*cols*3)/10;*/

  /* Mvec and Pvec arrays */
  bf_mvec_out=buffer_new(6,max_label,1);
  mvec_out=(momentvector *)bf_mvec_out->data;
  bf_pvec_out=buffer_new(3,max_label,1);
  pvec_out=bf_pvec_out->data;
  bf_mvec=buffer_new(6,max_label,1);
  mvec=(momentvector *)bf_mvec->data;
  bf_pvec=buffer_new(3,max_label,1);
  pvec=bf_pvec->data;
  
/* Other arrays */
  bf_aold=buffer_new(1,max_label,1);   /* Old areas */
  aold=bf_aold->data;
  bf_dalist=buffer_new(3,max_label,1);     /* rows: dist0,area0,int:(ts0,oind) */
  dalist=(da_tuple *)bf_dalist->data;
  bf_arate=buffer_new(6,max_label,1);      /* rows: (arate,amin,d0,dn,t0,tn) */
  arate=(arearate *)bf_arate->data;
  bf_bbox=ibuffer_new(4,max_label,1);      /* rows: xmin,xmax,ymin,ymax */
  bbox=(boundingbox *)bf_bbox->data;
  bf_nodes=ibuffer_new(2,max_label,1);     /* rows: prev,next */
  nodes=(node *)bf_nodes->data;
  active_list=NULL;
/*   active_end=NULL; */

  free_label=1; /* Running label number */
  boutcnt=0;    /* Number of blobs stored */

  d_pivot=0.0;  /* Current distance */
  timestep=0;   /* Count evolution steps */

  while(edges_left) {
    /* Move small valued edges to front */
    cnedge=0;
    while(cnedge==0) {
      if(tslist_flag) {
/* 	d_pivot=sqrt(thres[timestep]); */
	d_pivot=thres[timestep];
      } else {
	d_pivot+=res;
      }
      cnedge=edgelist_partition(&elist[first_edge],edges_left,d_pivot);
/* 				edges_left,SQR(d_pivot)); */
      timestep++;
    }
    /* Go through all edges in this time step */
    for(edgeno=first_edge;edgeno<first_edge+cnedge;edgeno++) {
      cx=elist[edgeno].xpos-1;  /* Remove matlab offset */
      cy=elist[edgeno].ypos-1;
      ind1=cy+rows*cx;
      switch((int)elist[edgeno].dirfl) {
      case 0:	/* Dx */
	ind2=ind1+rows;
	break;
      case 1:   /* Dy */	
	ind2=ind1+1;
	break;
      case 2:   /* Dxy */	
	ind2=ind1+1+rows;
	break;
      default:/* Dxy2 */	
	ind2=ind1+1-rows;
	break;
      }

      /* Read labels */
      l1=limg[ind1];
      l2=limg[ind2];
      if(l1==l2) {
	
      if(l1==0) {
	  /* New region */
	  cl=free_label++;    /* New label */

	  /* Insert at end of list */
/* 	  nodes[cl].next=NULL; */
/* 	  nodes[cl].prev=active_end; */
/* 	  if(active_end != NULL) active_end->next = &nodes[cl]; */
/* 	  active_end = &nodes[cl]; */
/* 	  if(active_list == NULL) active_list= &nodes[cl]; */

	  /* Insert at start of list */
	  nodes[cl].next = active_list;
	  nodes[cl].prev = NULL;
	  if(active_list != NULL) active_list->prev = &nodes[cl];
	  active_list = &nodes[cl];

	  if(free_label>max_label) {free_label=max_label; printf("OOPS!\n");}
	  /* Remember at what threshold we appeared */
	  dalist[cl].dist = d_pivot;
	  dalist[cl].area = 2;
	  dalist[cl].ts   = timestep;
	  dalist[cl].oind = 0;     /* Region not stored yet */
	  
      /* Set bounding box, mvec and pvec */
	  pvec_k = &pvec[cl*3];
	  ccx=cx+1;ccy=cy+1;   /* Add offset for moment computation */

      switch((int)elist[edgeno].dirfl) {
	  case 0:  /* Dx */
	    bbox[cl].xmin=cx;
	    bbox[cl].xmax=cx+1;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+(ccx+1);
	    mvec[cl].my=ccy+ccy;
	    mvec[cl].mx2=2*ccx*(ccx+1)+1;  /* cx*cx+(cx+1)*(cx+1) */
	    mvec[cl].mxy=2*ccx*ccy+ccy;    /* cx*cy+(cx+1)*cy */
	    mvec[cl].my2=2*ccy*ccy;
	    break;
	  case 1:  /* Dy */
	    bbox[cl].xmin=cx;
	    bbox[cl].xmax=cx;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy+1;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+ccx;
	    mvec[cl].my=ccy+(ccy+1);
	    mvec[cl].mx2=2*ccx*ccx; 
	    mvec[cl].mxy=2*ccx*ccy+ccx;     /* cx*cy+cx*(cy+1) */
	    mvec[cl].my2=2*ccy*(ccy+1)+1;   /* cy*cy+(cy+1)*(cy+1) */
	    break;
	  case 2:  /* Dxy */
	    bbox[cl].xmin=cx;
	    bbox[cl].xmax=cx+1;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy+1;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+(ccx+1);
	    mvec[cl].my=ccy+(ccy+1);
	    mvec[cl].mx2=2*ccx*(ccx+1)+1;    /* cx*cx+(cx+1)*(cx+1) */
	    mvec[cl].mxy=2*ccx*ccy+ccx+ccy+1;/* cx*cy+(cx+1)*(cy+1) */
	    mvec[cl].my2=2*ccy*(ccy+1)+1;    /* cy*cy+(cy+1)*(cy+1) */
	    break;
	  default: /* Dxy2 */
	    bbox[cl].xmin=cx-1;
	    bbox[cl].xmax=cx;
	    bbox[cl].ymin=cy;
	    bbox[cl].ymax=cy+1;
	    mvec[cl].area=2;
	    mvec[cl].mx=ccx+(ccx-1);
	    mvec[cl].my=ccy+(ccy+1);
	    mvec[cl].mx2=2*ccx*(ccx-1)+1;    /* (cx-1)*(cx-1)+cx*cx */
	    mvec[cl].mxy=2*ccx*ccy+ccx-ccy-1;/* cx*cy+(cx-1)*(cy+1) */
	    mvec[cl].my2=2*ccy*(ccy+1)+1;    /* cy*cy+(cy+1)*(cy+1) */
	  }
	  /* Initialize pvec */
	  pvec_k[0]=image[ind1]+image[ind2];
	  if(ndim>1) {
	    pvec_k[1]=image[ind1+rows*cols]+image[ind2+rows*cols];
	    pvec_k[2]=image[ind1+rows*cols*2]+image[ind2+rows*cols*2];
	  }
	  /* Set labels */
	  limg[ind1]=cl;
	  limg[ind2]=cl;      
	}
      } else {
	
     if(((l1==0)||(l2==0))&&(l1!=-1)&&(l2!=-1)) {
	  
      /* Append one pixel to a region */
	  if(l1>l2) {cl=l1;} else {cl=l2;} /* cl=max(l1,l2), non-zero label */
	  ccx=cx+1;ccy=cy+1;   /* Add offset for moment computation */
	  /* Update bounding box */
	  switch((int)elist[edgeno].dirfl) {
	  case 0:    /* Dx */
	    if(bbox[cl].xmin>cx)   bbox[cl].xmin=cx;
	    if(bbox[cl].xmax<cx+1) bbox[cl].xmax=cx+1;
	    if(cl==l1) ccx++;   /* Right one is new (cx+1,cy) */
	    break;
	  case 1:    /* Dy */
	    if(bbox[cl].ymin>cy)   bbox[cl].ymin=cy;
	    if(bbox[cl].ymax<cy+1) bbox[cl].ymax=cy+1;
	    if(cl==l1) ccy++;   /* Lower one is new (cx,cy+1) */
	    break;
	  case 2:    /* Dxy */
	    if(bbox[cl].ymin>cy)   bbox[cl].ymin=cy;
	    if(bbox[cl].ymax<cy+1) bbox[cl].ymax=cy+1;
	    if(bbox[cl].xmin>cx)   bbox[cl].xmin=cx;
	    if(bbox[cl].xmax<cx+1) bbox[cl].xmax=cx+1;
	    if(cl==l1) {ccx++;ccy++;}  /* Lower one is new (cx+1,cy+1) */
	    break;
	  default:   /* Dxy2 */
	    if(bbox[cl].ymin>cy)   bbox[cl].ymin=cy;
	    if(bbox[cl].ymax<cy+1) bbox[cl].ymax=cy+1;
	    if(bbox[cl].xmin>cx-1) bbox[cl].xmin=cx-1;
	    if(bbox[cl].xmax<cx)   bbox[cl].xmax=cx;
	    if(cl==l1) {ccx--;ccy++;}  /* Lower one is new (cx-1,cy+1) */
	  }
	  /* Update mvec */
	  mvec[cl].area += 1;
	  mvec[cl].mx   += ccx;
	  mvec[cl].my   += ccy;
	  mvec[cl].mx2  += ccx*ccx;
	  mvec[cl].mxy  += ccx*ccy;
	  mvec[cl].my2  += ccy*ccy;
	  
	  ind1=ccy-1+(ccx-1)*rows;   /* Index of new pixel */
	  /* Update pvec */
	  pvec_k = &pvec[cl*3];
	  pvec_k[0] += image[ind1];
	  if(ndim>1) {
	    pvec_k[1] += image[ind1+rows*cols];
	    pvec_k[2] += image[ind1+rows*cols*2];
	  }
	  /* Set label */
	  limg[ind1] = cl;
	} else {
	 
      if ( (l1!=-1) && (l2!=-1) )   
      {

      /* Merge regions */
	  /* Largest region keeps label */
 	  if(mvec[l1].area>mvec[l2].area) {cl=l1;lh=l2;} else {cl=l2;lh=l1;} 
	  /* Change all lh to cl */
	  pvec_k = &pvec[lh*3];
	  xmin=bbox[lh].xmin;
	  xmax=bbox[lh].xmax;
	  ymin=bbox[lh].ymin;
	  ymax=bbox[lh].ymax;
	  for(xx=xmin;xx<=xmax;xx++) {
	    for(yy=ymin;yy<=ymax;yy++) {
	      if(limg[yy+xx*rows]==lh) limg[yy+xx*rows]=cl;
	    }
	  }
	  /* Merge bounding boxes */
	  if(bbox[cl].xmin>bbox[lh].xmin) bbox[cl].xmin=bbox[lh].xmin;
	  if(bbox[cl].xmax<bbox[lh].xmax) bbox[cl].xmax=bbox[lh].xmax;
	  if(bbox[cl].ymin>bbox[lh].ymin) bbox[cl].ymin=bbox[lh].ymin;
	  if(bbox[cl].ymax<bbox[lh].ymax) bbox[cl].ymax=bbox[lh].ymax;
	  /* Merge mvecs */
	  mvec[cl].area += mvec[lh].area;
	  mvec[cl].mx   += mvec[lh].mx;
	  mvec[cl].my   += mvec[lh].my;
	  mvec[cl].mx2  += mvec[lh].mx2;
	  mvec[cl].mxy  += mvec[lh].mxy;
	  mvec[cl].my2  += mvec[lh].my2;
	  /* Merge pvecs */
	  pvec_l = &pvec[cl*3];
	  pvec_l[0] += pvec_k[0];
	  pvec_l[1] += pvec_k[1];
	  pvec_l[2] += pvec_k[2];
	  /* Indicate as merged */
	  mvec[lh].area=0;

	  /* Unlink node from list */
	  if(nodes[lh].prev != NULL) 
	    nodes[lh].prev->next = nodes[lh].next;
	  else
	    active_list = nodes[lh].next; /* First in list removed */
	  if(nodes[lh].next != NULL) 
	    nodes[lh].next->prev = nodes[lh].prev;
/* 	  else */
/* 	    active_end = nodes[lh].prev; /\* Last in list removed *\/ */
	}}
      }
    }
    
   /*
     * New step in evolution, check all regions
     */
/*     cl=limg[285+215*rows]; */
/*     new_arate=(mvec[cl].area-dalist[cl].area)/(d_pivot-dalist[cl].dist); */
/*     printf("%d: dalist=[%g %g %d] nfl=%d\tarate=%g \ta=%g ao=%g cfl=%d\n",timestep,dalist[cl].area,dalist[cl].dist,dalist[cl].ts,mvec[cl].area>aold[cl]*size_inc,new_arate,mvec[cl].area,aold[cl],mvec[cl].area<dalist[cl].area*2.0); */
/*     for(k=0;k<top_label;k++) { */
    cnode = active_list;
    while(cnode != NULL) {
      /* Get index */
      k = cnode-nodes;
      cnode = cnode->next; /* Next node to visit */
      /* Require min area */
      if(mvec[k].area>amin) {
	/* Require small area increase since start */
	if(mvec[k].area<dalist[k].area*2.0) {
/*   	  new_arate=mvec[k].area/dalist[k].area/(d_pivot-dalist[k].dist);  */
	  if(timestep-dalist[k].ts>2) {
	    new_arate=(mvec[k].area-dalist[k].area)/(d_pivot-dalist[k].dist); 
	  } else {
	    new_arate=1e14;
	  }
	  /* Find out where to store region */
	  out_ind=dalist[k].oind;
	  if(!out_ind) {
	    /* New out_ind */
	    out_ind=boutcnt++;
	    dalist[k].oind=out_ind;
	    /* Remember start posn */
	    arate[out_ind].d0 = dalist[k].dist; /* d_0 */
	    arate[out_ind].t0 = dalist[k].ts;   /* t_0 */
	    old_arate=1e15;
	  } else {
	    old_arate=arate[out_ind].arate;
	  }
	  if(new_arate<old_arate) { 
	    /* Update blob description */
	    mvec_out[out_ind]=mvec[k];
	    memcpy(&pvec_out[out_ind*3],&pvec[k*3],3*sizeof(fpnum));
	    /* Update opt point */
	    arate[out_ind].arate = new_arate;
	    arate[out_ind].amin  = mvec[k].area;
 	  }
	  /* Update end point */
	  arate[out_ind].dn    = d_pivot;
 	  arate[out_ind].tn    = timestep; 
	} else {
	  if(dalist[k].oind) {
	    if(mvec[k].area<aold[k]*size_inc) {
	      mvec_out[dalist[k].oind].area=0;         /* Remove region */
	    }
	  }
	}
	if(boutcnt>max_label) {
	  printf("OOPS!");
	  boutcnt=max_label;
	}
      }
      /* Check if large size increase */
      if(mvec[k].area>aold[k]*size_inc) {
	/* NOTE: We cannot set end point here, this will miss mergers */
	/* NOTE: Regions start with area=2, then get correct area here */
	dalist[k].dist = d_pivot;      /* Reset start distance */
	dalist[k].area = mvec[k].area; /* Reset start area */
	dalist[k].ts   = timestep;     /* Reset start timestep */
	dalist[k].oind = 0;            /* Terminate */
      }
      aold[k]=mvec[k].area; /* Set aold */
    }
    /* End of while loop */
    first_edge += cnedge;
    edges_left -= cnedge;
      
/*      printf("%d ",cnedge);  */
/*     if(cnedge>0) { */
/*       update_edgelist(bf_image,&elist[first_edge], */
/* 		      edges_left,bf_limg,bf_pvec,bf_mvec); */
/*     } */
  }
  ibuffer_free(bf_bbox);
  ibuffer_free(bf_limg);
  ibuffer_free(bf_nodes);
/*   printf("boutcnt=%d\n",boutcnt); */
/*   printf("timestep=%d d_pivot=%g\n",timestep,d_pivot); */
  buffer_free(bf_mvec);
  buffer_free(bf_aold);
  buffer_free(bf_dalist);
  buffer_free(bf_pvec);
  bfp_mvec[0]=bf_mvec_out;
  bfp_pvec[0]=bf_pvec_out;
  bfp_arate[0]=bf_arate;
}

