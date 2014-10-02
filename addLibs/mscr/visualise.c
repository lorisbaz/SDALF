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

#include "visualise.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

/*
** Paint an image buffer in a given colour
**
** bf       image buffer
** bf_pvec  desired colour
**
*/
void buffer_paint(buffer *bf,buffer *bf_pvec) {
  int rows,cols,ndim;
  int k,l,m;
  fpnum *img,ccol;

  rows=bf->rows;
  cols=bf->cols;
  ndim=bf->ndim;
  img =bf->data;

  for(m=0;m<ndim;m++) {
    ccol=bf_pvec->data[m];
    for(l=0;l<cols;l++) {
      for(k=0;k<rows;k++) {
	img[k+l*rows+m*rows*cols]=ccol;
      }
    }
  }
}
/*
** Paint an image byte buffer in a given colour
**
** bf       image byte buffer
** bf_pvec  desired colour
**
*/
void bbuffer_paint(bbuffer *bf,bbuffer *bf_pvec) {
  int rows,cols,ndim;
  int k,l,m;
  unsigned char *img,ccol;

  rows=bf->rows;
  cols=bf->cols;
  ndim=bf->ndim;
  img =bf->data;

  for(m=0;m<ndim;m++) {
    ccol=bf_pvec->data[m];
    for(l=0;l<cols;l++) {
      for(k=0;k<rows;k++) {
	img[k+l*rows+m*rows*cols]=ccol;
      }
    }
  }
}
/*
 * Closed form eigenvalue decomposition
 *
 * I  Positive definite input matrix of form
 *       I[0] I[1]
 *       I[2] I[3]   where I[1]=I[2]
 * D  Output eigenvalues  D[0]>=D[1]
 * E  Eigenvector matrix of form
 *       E[0] E[2]
 *       E[1] E[3]
 */
void eigendec(fpnum *I,fpnum *D,fpnum *E) {
fpnum a,b,enorm;
 
 /* Find eigenvalues */ 
 a=I[0]+I[3];                         /* trace */
 b=sqrt(SQR(I[0]-I[3])+4*SQR(I[1]));  
 D[0]=(a+b)/2;
 D[1]=(a-b)/2;
 
 /* Find eigenvector 1 */ 
 E[0]=I[3]+I[1]-D[0];
 E[1]=D[0]-I[0]-I[1];
 enorm=sqrt(SQR(E[0])+SQR(E[1]));
 E[0]=E[0]/enorm;
 E[1]=E[1]/enorm;

 /* Find eigenvector 2 */ 
 E[2]=I[3]+I[1]-D[1];
 E[3]=D[1]-I[0]-I[1];
 enorm=sqrt(SQR(E[2])+SQR(E[3]));
 E[2]=E[2]/enorm;
 E[3]=E[3]/enorm;
}
/*
** Draw a filled ellipse in an RGB image.
**
** bf_img  Image buffer to paint in
** mvec    Moment list (6x1)
** pvec    Property (colour) list (3x1)
**
*/
void draw_ellipse3(buffer *bf_img,fpnum *mvec,fpnum *pvec)
{
int x,y,k;
int rows,cols,pixels,ndim;
int xl,xh,yl,yh;
int cind;
fpnum xmax,ymax,xd,yd,qf,qf0,xd2a12;
fpnum cogx,cogy;
fpnum cov[]={0,0,0,0};
fpnum D[]={0,0};
fpnum E[]={0,0,0,0};
fpnum a11,a12,a22;
fpnum *img;

  rows=bf_img->rows;
  cols=bf_img->cols;
  ndim=bf_img->ndim;
  pixels=rows*cols;
  img=bf_img->data;

  /* Extract centroid and covariance matrix */
  cogx=mvec[1];
  cogy=mvec[2];
  cov[0]=mvec[3];
  cov[1]=mvec[4];
  cov[2]=mvec[4];
  cov[3]=mvec[5];
  /* Decompose covariance matrix */
  eigendec(cov,D,E);

  xmax=2.0*sqrt(D[0]*SQR(E[0])+D[1]*SQR(E[2]));
  /* Interval to check, with matlab offset removed */
  xl=(int)ceil(cogx-xmax-1.0);if(xl<0) xl=0;
  xh=(int)floor(cogx+xmax-1.0);if(xh>cols-1) xh=cols-1;
  ymax=2.0*sqrt(D[0]*SQR(E[1])+D[1]*SQR(E[3]));
  /* Interval to check, with matlab offset removed */
  yl=(int)ceil(cogy-ymax-1.0);if(yl<0) yl=0;
  yh=(int)floor(cogy+ymax-1.0);if(yh>rows-1) yh=rows-1;

  /* A= .25*E*D^-1*E' */
  a11=.25*(E[0]*E[0]/D[0]+E[2]*E[2]/D[1]);
  a12=.25*(E[1]*E[0]/D[0]+E[2]*E[3]/D[1]);
  /* a21=.25*(E[0]*E[1]/D[0]+E[2]*E[3]/D[1]);*/
  a22=.25*(E[1]*E[1]/D[0]+E[3]*E[3]/D[1]);

  /* (x-m)'*A*(x-m)<1 */
  for(x=xl;x<=xh;x++) {
    xd=(fpnum)x-cogx+1.0;   /* Add matlab offset again */
    qf0=xd*xd*a11;
    xd2a12=xd*2*a12;
    for(y=yl;y<=yh;y++) {
      yd=(fpnum)y-cogy+1.0; /* Add matlab offset again */
      qf=qf0;
      qf+=yd*(xd2a12+yd*a22);
      if(qf<1.0) {
	cind=x*rows+y;
	for(k=0;k<ndim;k++) {
	  img[cind+k*pixels]=pvec[k];
	}
      }
    }
  }
}
/*
** Draw a filled ellipse in an uint8 RGB image.
**
** bf_img  Image buffer to paint in
** mvec    Moment list (6x1)
** pvec    Property (colour) list (3x1)
**
*/
void bdraw_ellipse3(bbuffer *bf_img,fpnum *mvec,unsigned char *pvec)
{
int x,y,k;
int rows,cols,pixels,ndim;
int xl,xh,yl,yh;
int cind;
fpnum xmax,ymax,xd,yd,qf,qf0,xd2a12;
fpnum cogx,cogy;
fpnum cov[]={0,0,0,0};
fpnum D[]={0,0};
fpnum E[]={0,0,0,0};
fpnum a11,a12,a22;
unsigned char *img;

  rows=bf_img->rows;
  cols=bf_img->cols;
  ndim=bf_img->ndim;
  pixels=rows*cols;
  img=bf_img->data;

  /* Extract centroid and covariance matrix */
  cogx=mvec[1];
  cogy=mvec[2];
  cov[0]=mvec[3];
  cov[1]=mvec[4];
  cov[2]=mvec[4];
  cov[3]=mvec[5];
  /* Decompose covariance matrix */
  eigendec(cov,D,E);

  xmax=2.0*sqrt(D[0]*SQR(E[0])+D[1]*SQR(E[2]));
  /* Interval to check, with matlab offset removed */
  xl=(int)ceil(cogx-xmax-1.0);if(xl<0) xl=0;
  xh=(int)floor(cogx+xmax-1.0);if(xh>cols-1) xh=cols-1;
  ymax=2.0*sqrt(D[0]*SQR(E[1])+D[1]*SQR(E[3]));
  /* Interval to check, with matlab offset removed */
  yl=(int)ceil(cogy-ymax-1.0);if(yl<0) yl=0;
  yh=(int)floor(cogy+ymax-1.0);if(yh>rows-1) yh=rows-1;

  /* A= .25*E*D^-1*E' */
  a11=.25*(E[0]*E[0]/D[0]+E[2]*E[2]/D[1]);
  a12=.25*(E[1]*E[0]/D[0]+E[2]*E[3]/D[1]);
  /* a21=.25*(E[0]*E[1]/D[0]+E[2]*E[3]/D[1]);*/
  a22=.25*(E[1]*E[1]/D[0]+E[3]*E[3]/D[1]);

  /* (x-m)'*A*(x-m)<1 */
  for(x=xl;x<=xh;x++) {
    xd=(fpnum)x-cogx+1.0;   /* Add matlab offset again */
    qf0=xd*xd*a11;
    xd2a12=xd*2*a12;
    for(y=yl;y<=yh;y++) {
      yd=(fpnum)y-cogy+1.0; /* Add matlab offset again */
      qf=qf0;
      qf+=yd*(xd2a12+yd*a22);
      if(qf<1.0) {
	cind=x*rows+y;
	for(k=0;k<ndim;k++) {
	  img[cind+k*pixels]=pvec[k];
	}
      }
    }
  }
}
/*
** Sort blobs according to size, and visualise
** as ellipses in an RGB image.
**
** bf_img  Image buffer to paint in
** bf_mvec Moment list (6xN)
** bf_pvec Property (colour) list (3xN)
**
*/
void draw_ellipses(buffer *bf_img,buffer *bf_mvec,buffer *bf_pvec) {
int k,l,m;
int nblobs; 
int *indl;
fpnum *mvec_data,*pvec_data;

 nblobs=bf_mvec->cols; 
 mvec_data=bf_mvec->data;
 pvec_data=bf_pvec->data;

 /* Allocate index list */  
 indl=(int *)calloc(nblobs,sizeof(int)); 

 /*
 ** Fill in index list using
 ** insertion sort of blob areas
 */

 indl[0]=0;  /* Insert first index */
 for(k=1;k<nblobs;k++) {
   /* Find insertion point */
   l=0;
   while((mvec_data[indl[l]*6]>mvec_data[k*6])&&(l<k)) l++;
   /* Move old indices up */
   for(m=k;m>l;m--) indl[m]=indl[m-1];
   /* insert new index */
   indl[l]=k;
 } 
 for(k=0;k<nblobs;k++) {
   draw_ellipse3(bf_img,&mvec_data[indl[k]*6],&pvec_data[indl[k]*3]);
 }
 free(indl);
}
/*
** Sort blobs according to size, and visualise
** as ellipses in an uint8 RGB image.
**
** bf_img  Image buffer to paint in
** bf_mvec Moment list (6xN)
** bf_pvec Property (colour) list (3xN)
**
*/
void bdraw_ellipses(bbuffer *bf_img,buffer *bf_mvec,bbuffer *bf_pvec) {
int k,l,m;
int nblobs;
int *indl;
fpnum *mvec_data;
unsigned char *pvec_data;

 nblobs=bf_mvec->cols; 
 mvec_data=bf_mvec->data;
 pvec_data=bf_pvec->data;

 /* Allocate index list */  
 indl=(int *)calloc(nblobs,sizeof(int)); 

 /*
 ** Fill in index list using
 ** insertion sort of blob areas
 */

 indl[0]=0;  /* Insert first index */
 for(k=1;k<nblobs;k++) {
   /* Find insertion point */
   l=0;
   while((mvec_data[indl[l]*6]>mvec_data[k*6])&&(l<k)) l++;
   /* Move old indices up */
   for(m=k;m>l;m--) indl[m]=indl[m-1];
   /* insert new index */
   indl[l]=k;
 } 
 for(k=0;k<nblobs;k++) {
   bdraw_ellipse3(bf_img,&mvec_data[indl[k]*6],&pvec_data[indl[k]*3]);
 }
 free(indl);
}
/*
** Paint regions with their average colours
** in an RGB image.
**
** bf_img     Image buffer to paint in
** bf_labelim Region label image
** bf_pvec    Property (colour) list (3xN)
**
 */
void draw_regions(buffer *bf_img,ibuffer *bf_labelim,buffer *bf_pvec) {
int x,y,k;
int cind,cl;
int nblobs,rows,cols,ndim,pixels; 
fpnum *img,*pvec_data;
int *labelim;

 rows=bf_img->rows;
 cols=bf_img->cols;
 ndim=bf_img->ndim;
 pixels=rows*cols;
 img=bf_img->data;

 labelim=bf_labelim->data;

 nblobs=bf_pvec->cols;    /* Nblobs after merging */
 pvec_data=bf_pvec->data;

 /* Loop over image and paint */
 for(x=0;x<cols;x++) {
   for(y=0;y<rows;y++) {
     cind=y+x*rows;
     cl=labelim[cind];
     if((cl>0)&&(cl<=nblobs)) {
       for(k=0;k<ndim;k++) {
	 img[cind+k*pixels]=pvec_data[(cl-1)*ndim+k];
       }
     }
   }
 }
}
/*
** Convert a greyscale pvec to RGB by copying the
** grey-value 3x.
**
** bf_pvec    Input pvec (1xN)
** bl_pvecn   Pointer to output pvec (3xN)
**
 */
void pvec_to_rgb(buffer *bf_pvec,buffer **bl_pvecn) {
int rows,cols,k;
buffer *bf_pvecn;
fpnum *data0,*data1;
 rows=bf_pvec->rows;
 cols=bf_pvec->cols;
 data0=bf_pvec->data;
 if(rows!=1) {
   printf("pvec_to_rgb:Error bf_pvec is not 1xN\n");
   return;
 }
 bf_pvecn=buffer_new(3,cols,1);
 data1=bf_pvecn->data;
 for(k=0;k<cols;k++) {
   data1[k*3]=data1[k*3+1]=data1[k*3+2]=data0[k];
 }
 bl_pvecn[0]=bf_pvecn;
}
/*
 *  Change back to RGB space
 */
void pvec_colourspace(buffer *bf_pvec) {
  int rows,cols,k;
  fpnum *data,v1,v2,v3,r,g,b;
 rows=bf_pvec->rows;
 cols=bf_pvec->cols;
 data=bf_pvec->data;
 for(k=0;k<cols;k++) {
   v1=data[k*3];      /* (r+g+b)/3 */
   v2=data[k*3+1]; /* r-g */
   v3=data[k*3+2]; /* (r+g)/2-b */
   /* b= v1-v3*2/3 */
   b=v1-v3*2.0/3.0;
   /* r=v3+b+v2/2 */
   r= v3+b+v2/2.0;
   /* g= 3*v1-r-b*/
   g=3*v1-r-b;
   data[k*3]=r;
   data[k*3+1]=g;
   data[k*3+2]=b;
 }
}
/*
** Convert an RGB fpnum pvec to uint8.
**
** bf_pvec    Input fpnum pvec (3xN)
** bl_pvec2  Output uint8 pvec (3xN)
**
 */
void pvec_to_uint8(buffer *bf_pvec,bbuffer **bl_pvec2) {
  int rows,cols,k,m;
bbuffer *bf_pvec2;
fpnum *data0;
unsigned char *data1;
 rows=bf_pvec->rows;
 cols=bf_pvec->cols;
 data0=bf_pvec->data;
 bf_pvec2=bbuffer_new(rows,cols,1);
 data1=bf_pvec2->data;
 for(k=0;k<cols;k++) {
   for(m=0;m<rows;m++) {
     data1[m+k*rows]=(unsigned char)(data0[m+k*rows]*255.0);
   }
 }
 bl_pvec2[0]=bf_pvec2;
}
/*
** Destructively set area field in moment vector
** to the area of the approximating ellipse.
**
** bf_mvec    Input mvec (6xN)
**
 */
void mvec_set_area(buffer *bf_mvec) {
int cols,k,cind;
fpnum *mvec;
 cols = bf_mvec->cols;
 mvec = (fpnum *)bf_mvec->data;
 for(k=0;k<cols;k++) {
   cind=k*6;
   mvec[cind] = 4.0*PI*sqrt(mvec[cind+3]*mvec[cind+5]-mvec[cind+4]*mvec[cind+4]);
 }
}
/*
 * Compute ratios of mask pixels inside ellipses
 *
 * bf_mask    Binary mask
 * bf_mvec    Input mvec (6xN)
 * bf_ratios  Resultant ratios
 *
 */
void inside_mask_ratios(ibuffer *bf_mask,buffer *bf_mvec,buffer *bf_ratios) {
int nblobs,k,area_e,area_em,rows,cols;
fpnum *ratios,*mvec,*cmvec;
int *mask;
int x,y,xl,xh,yl,yh;
fpnum xmax,ymax,xd,yd,qf,qf0,xd2a12;
fpnum cogx,cogy;
fpnum cov[]={0,0,0,0};
fpnum D[]={0,0};
fpnum E[]={0,0,0,0};
fpnum a11,a12,a22;

  nblobs = bf_mvec->cols; 
  mvec   = bf_mvec->data;
  ratios = bf_ratios->data;
  mask   = bf_mask->data;
  rows   = bf_mask->rows;
  cols   = bf_mask->cols;

  for(k=0;k<nblobs;k++) {
    /* Reset area counts */
    area_em=area_e=0;
    
    /* Extract centroid and covariance matrix */
    cmvec  = &mvec[k*6];
    cogx   = cmvec[1];
    cogy   = cmvec[2];
    cov[0] = cmvec[3];
    cov[1] = cmvec[4];
    cov[2] = cmvec[4];
    cov[3] = cmvec[5];
    /* Decompose covariance matrix */
    eigendec(cov,D,E);
    
    xmax=2.0*sqrt(D[0]*SQR(E[0])+D[1]*SQR(E[2]));
    /* Interval to check, with matlab offset removed */
    xl=(int)ceil(cogx-xmax-1.0);if(xl<0) xl=0;
    xh=(int)floor(cogx+xmax-1.0);if(xh>cols-1) xh=cols-1;
    ymax=2.0*sqrt(D[0]*SQR(E[1])+D[1]*SQR(E[3]));
    /* Interval to check, with matlab offset removed */
    yl=(int)ceil(cogy-ymax-1.0);if(yl<0) yl=0;
    yh=(int)floor(cogy+ymax-1.0);if(yh>rows-1) yh=rows-1;
    
    /* A= .25*E*D^-1*E' */
    a11=.25*(E[0]*E[0]/D[0]+E[2]*E[2]/D[1]);
    a12=.25*(E[1]*E[0]/D[0]+E[2]*E[3]/D[1]);
    /* a21=.25*(E[0]*E[1]/D[0]+E[2]*E[3]/D[1]);*/
    a22=.25*(E[1]*E[1]/D[0]+E[3]*E[3]/D[1]);
    
    /* (x-m)'*A*(x-m)<1 */
    for(x=xl;x<=xh;x++) {
      xd=(fpnum)x-cogx+1.0;   /* Add matlab offset again */
      qf0=xd*xd*a11;
      xd2a12=xd*2*a12;
      for(y=yl;y<=yh;y++) {
	yd=(fpnum)y-cogy+1.0; /* Add matlab offset again */
	qf=qf0;
	qf+=yd*(xd2a12+yd*a22);
	if(qf<1.0) {
	  /* We're now inside ellipse */
	  area_e++;
	  if(mask[x*rows+y]>0) area_em++;
	}
      }
    }
    
    /* Store ration of areas */
    ratios[k]=(fpnum)area_em/area_e;
  }
}
