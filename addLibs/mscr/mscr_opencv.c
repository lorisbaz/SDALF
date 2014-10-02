/*
**  Demo of MSCR using OpenCV GUI (highgui)
**
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

#include <stdlib.h>
#include <stdio.h>

#include "cv.h"
#include "highgui.h"

#include "measure_time.h"
#include "image_buffer.h"
#include "msr_util.h"
#include "visualise.h"

/* Default settings */
#define DEF_MARGIN 0.0015
#define DEF_TIMESTEPS 200
#define DEF_N8FLAG 0
#define DEF_NORMFL 1
#define DEF_BLURFL 0
#define DEF_LOWRESFL 1
#define DEF_TIMESIZE 10
#define DEF_CSPACE 0
/*
 *  Global variables
 */

CvCapture *cam;        /* Highgui camera handle */
double tb_margin;
int tb_timesteps;
int kb_n8flag=DEF_N8FLAG;     /* 1 for 8 neighbours instead of 4 */
int kb_normfl=DEF_NORMFL;     /* 1 for Chi2 edges instead of Euclidean */
int kb_blurfl=DEF_BLURFL;     /* 1 for image blur instead of edge blur */
int kb_lowresfl=DEF_LOWRESFL; /* Toggles 320x240 or 640x480 */
int kb_cspace=DEF_CSPACE;     /* Toggles Y,Cr,Cb and RGB */

#define TB1_MAX 2000
#define TB1_SCALE 0.05
#define TB2_MAX 1000

/*
 *  Callback functions
 */
void callback_trackbar1(int cval) {
  tb_margin = (double)cval/TB1_MAX*TB1_SCALE;
  printf("min_margin=%g\n",tb_margin);
};

void callback_trackbar2(int cval) {
  tb_timesteps=cval;
  printf("timesteps =%d\n",tb_timesteps);
};

void callback_mouse(int event, int x, int y, int flags, void* param) {
  if(event==CV_EVENT_LBUTTONDBLCLK) {
    cvDestroyAllWindows();
    cvReleaseCapture(&cam);
    exit(0);
  }
}

/*
 *  Copy and flip contents of image1 to image2
 *  If image2 is bigger than image1, only the upper left corner is used
 */
void flip_image(IplImage *image1,IplImage *image2) {
unsigned char *image1_data, *image2_data;
 int width,height,ndim,x,y,m;
 int width2,height2;

 image1_data = (unsigned char *)image1->imageData;
 image2_data = (unsigned char *)image2->imageData;
 height = image1->height;
 width = image1->width;
 height2 = image2->height;
 width2 = image2->width;
 ndim = image1->nChannels;
 for(y=0;y<height;y++) {
   for(x=0;x<width;x++) {
     for(m=0;m<ndim;m++) {
       image2_data[m+(width-1-x)*ndim+y*width2*ndim]=
	 image1_data[m+x*ndim+y*width*ndim];
     }
   }
 }
}
/*
 *  Copy and contents of image1 to image2 and upsample
 */
void upsample_image(IplImage *image1,IplImage *image2) {
unsigned char *image1_data, *image2_data;
 int width,height,ndim,x,y,m;
 int width2,height2;
 int cind;
 unsigned char cval;

 image1_data = (unsigned char *)image1->imageData;
 image2_data = (unsigned char *)image2->imageData;
 height = image1->height;
 width = image1->width;
 height2 = image2->height;
 width2 = image2->width;
 ndim = image1->nChannels;
 for(y=0;y<height;y++) {
   for(x=0;x<width;x++) {
     for(m=0;m<ndim;m++) {
       cval = image1_data[m+x*ndim+y*width*ndim];
       cind = m+x*2*ndim+y*2*width2*ndim;
       image2_data[cind]= cval;
       image2_data[cind+ndim]= cval;
       image2_data[cind+width2*ndim]= cval;
       image2_data[cind+ndim+width2*ndim]= cval;
     }
   }
 }
}

/*
**  Copy an ipl image into a buffer struct, overwriting contents
**  optionally flipping the x-coordinate
*/
void buffer_iplcopy(buffer *bf_img,IplImage *ipl_img,int flipfl) {
  fpnum *img;
  unsigned char *iplimg_data;
  int width,height,ndim,x,y,m;
  height = ipl_img->height;
  width = ipl_img->width;
  ndim = ipl_img->nChannels;
  if(height != bf_img->rows) {
    printf("Error: image and buffer have different sizes!\n");
    printf("image size=%dx%dx%d, buffer size=%dx%dx%d\n",height,width,ndim,bf_img->rows,bf_img->cols,bf_img->ndim);
    printf("Maybe cvSetCaptureProperty is ignored? Try setting DEF_LOWRESFL 0 and recompile.\n");
    exit(1);
  }
  img = bf_img->data;
  iplimg_data = (unsigned char *)ipl_img->imageData;
  if(flipfl) {
    for(m=0;m<ndim;m++) {
      for(x=0;x<width;x++) {
	for(y=0;y<height;y++) {
	  img[y+x*height+m*width*height] = 
	    (fpnum)iplimg_data[m+(width-1-x)*ndim+y*width*ndim]/255.0;
	}
      }
    }
  } else {
    for(m=0;m<ndim;m++) {
      for(x=0;x<width;x++) {
	for(y=0;y<height;y++) {
	  img[y+x*height+m*width*height] = 
	    (fpnum)iplimg_data[m+x*ndim+y*width*ndim]/255.0;
	}
      }
    }
  }
}

/*
**  Copy an image buffer into an ipl image, overwriting contents
*/
void paste_image(IplImage *ipl_img,bbuffer *bf_img,int xoffset) {
  unsigned char *img;
  unsigned char  *iplimg_data;
  int width,height,ndim,x,y,m;
  int width2,height2;

 img = bf_img->data;
 height = bf_img->rows;
 width = bf_img->cols;
 ndim = bf_img->ndim;
 iplimg_data = (unsigned char *)ipl_img->imageData;
 height2 = ipl_img->height;
 width2 = ipl_img->width;
 for(m=0;m<ndim;m++) {
   for(x=0;x<width;x++) {
     for(y=0;y<height;y++) {
       iplimg_data[m+(x+xoffset)*ndim+y*width2*ndim]=
	 img[y+x*height+m*width*height];
     }
   }
 }
}
/*
**  Main function
*/
void detect_and_render_blobs(buffer *bf_image, bbuffer *bf_blobimg) {
buffer *bf_pvec,*bf_pvec2;
bbuffer *bf_pvec8;
ebuffer *bf_elist,*bf_thres=NULL;
buffer *bf_mvec,*bf_mvec2,*bf_arate;
buffer *bf_image2;
int validcnt,nofedges,rows,cols,ndim;
edgeval d_max;
double d_mean;
/* For draw_blobs */
//fpnum bkgr_green[] = {0.0,1.0,0.0}; /* Background colour */
unsigned char bkgr_green[] ={0,255,0}; /* Background colour */
bbuffer *bf_bkgr;

int min_size=60;          /* Default */
double ainc=1.01;         /* Default */
edgeval min_margin;       /* Value after conversion */
fpnum res=1e-4;           /* Default */
int timesteps=200;        /* Default */
int filter_size=3;        /* Default */
int n8flag=1;             /* Default */
int normfl=1;             /* Default */
int blurfl=0;             /* Default */
int verbosefl=0;          /* Default */

/* Set margin */
 min_margin  = tb_margin*EDGE_SCALE+EDGE_OFFSET;
 
 /* Set timesteps */
 timesteps = tb_timesteps;

 /* Set n8flag */
 n8flag=kb_n8flag;

 /* Set n8flag */
 normfl=kb_normfl;

 /* Set blurfl */
 blurfl=kb_blurfl;

/* Set edge list size */
 rows = bf_image->rows;
 cols = bf_image->cols;
 ndim = bf_image->ndim;
 if(n8flag) 
   nofedges=4*rows*cols-3*rows-3*cols+2;
 else
   nofedges=2*rows*cols-rows-cols;

 if(kb_cspace) {
   bf_image2=buffer_new(bf_image->rows,bf_image->cols,bf_image->ndim);
   image_colourspace(bf_image,bf_image2);
   bf_image=bf_image2;
 }

 bf_elist = ebuffer_new(4,nofedges,1);

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
   printf("bfz should be odd.\n");
 }

 /* Call thresholds interpolation */
 bf_thres=ebuffer_new(1,timesteps,1);
 d_mean=evolution_thresholds2(bf_elist,bf_thres,ndim);
 if(verbosefl) printf("d_mean=%g\n",d_mean);
 /* Try linear dependence on mean edge strength */
 /*      min_margin=d_mean*min_margin_sc; */

 /* Call computation function */
 edgelist_to_bloblist(&bf_mvec,&bf_pvec,&bf_arate,bf_image,bf_elist,bf_thres,min_size,ainc,res,verbosefl);
 center_moments(bf_mvec,bf_pvec);
 validcnt=bloblist_mark_invalid(bf_mvec,min_size,bf_arate,(fpnum)min_margin);
 validcnt=bloblist_shape_invalid(bf_mvec);
 if(verbosefl) printf("validcnt=%d\n",validcnt);
 
 /* Release timestep list */
 ebuffer_free(bf_thres);

 /* Allocate out arrays */
 bf_mvec2 = buffer_new(6,validcnt,1);
 bf_pvec2 = buffer_new(3,validcnt,1);
 
 bloblist_compact2(bf_mvec,bf_mvec2,bf_pvec,bf_pvec2);

 buffer_free(bf_mvec);   /* Release non-compacted mvec */
 buffer_free(bf_pvec);   /* Release non-compacted pvec */
 buffer_free(bf_arate);  /* Release non-compacted arate */

 /* Create an empty green image */
 bf_bkgr=bbuffer_new0(bkgr_green,3,1,1);  /* Use default background */
 bbuffer_paint(bf_blobimg,bf_bkgr);
   
 /* Set area to approximating ellipse area */
 mvec_set_area(bf_mvec2);  /* Blobs are drawn in descending area order */

 /* Convert pvec to uint8 */
 if(kb_cspace) pvec_colourspace(bf_pvec2);
 pvec_to_uint8(bf_pvec2,&bf_pvec8);

 /* Draw blobs on top of background */
 bdraw_ellipses(bf_blobimg,bf_mvec2,bf_pvec8);

 /* release memory */
 ebuffer_free(bf_elist);
 buffer_free(bf_mvec2);
 buffer_free(bf_pvec2);
 bbuffer_free(bf_pvec8);
 free(bf_bkgr);
 if(kb_cspace) buffer_free(bf_image);
}

int main( int argc, char** argv ) {
  int ikey=0;
  char key;
  IplImage *imCapt,*imDisp,*imDisp2;
  static int tb1_val,tb2_val;  /* Trackbar parameters */
  timewin *tw;
  buffer *bf_img;
  bbuffer *bf_blobs;
  int width,height;

  /* Set window title to a list of control keys */
  const char *window1 = "Controls: (n)eighbours, (d)istance, (b)lur type, (h)igh-res, (c)olourspace, (r)eset, (ESC) exit";
  const char *trackbar1 = "margin";
  const char *trackbar2 = "timesteps";
  
  int fno;
 
  if(kb_lowresfl) {
    width=320;
    height=240;
  } else {
    width=640;
    height=480;
  }

  cvInitSystem( argc,argv );

  /* Get an OpenCV camera handle */
  //cam = cvCreateCameraCapture(-1);
  cam = cvCreateCameraCapture(0);

  if(cam == NULL) {
    fprintf(stderr,"No camera found. Terminating.\n");
    exit(1);
  }

  /* Set size of image (appears to be ignored in Linux) */
  cvSetCaptureProperty(cam, CV_CAP_PROP_FRAME_WIDTH, width);
  cvSetCaptureProperty(cam, CV_CAP_PROP_FRAME_HEIGHT, height);

  /* Create a window with slider */
  cvNamedWindow(window1, CV_WINDOW_AUTOSIZE);
  cvSetMouseCallback(window1,callback_mouse, NULL );
  tb_margin = DEF_MARGIN; /* Default */
  tb1_val=tb_margin/TB1_SCALE*TB1_MAX;
  cvCreateTrackbar(trackbar1,window1,&tb1_val,TB1_MAX,callback_trackbar1);
  tb2_val = tb_timesteps = DEF_TIMESTEPS; /* Default */
  cvCreateTrackbar(trackbar2,window1,&tb2_val,TB2_MAX,callback_trackbar2);
  cvMoveWindow(window1, 100, 45);

  /* Allocate image buffers */
  if(kb_lowresfl)
    imDisp2 =  cvCreateImage( cvSize(width*4,height*2), IPL_DEPTH_8U, 3);
  else
    imDisp2 =  cvCreateImage( cvSize(width*2,height), IPL_DEPTH_8U, 3);
  imDisp =  cvCreateImage( cvSize(width*2,height), IPL_DEPTH_8U, 3);
  bf_img = buffer_new(height,width,3);
  bf_blobs = bbuffer_new(height,width,3);

  tw=timewin_new(DEF_TIMESIZE);
  fno=0;

  key=(char)cvWaitKey(500);
  
  while ( key !=27 ) {

    imCapt = cvQueryFrame(cam);

    buffer_iplcopy(bf_img,imCapt,1);

    /* Detect blobs */
    detect_and_render_blobs(bf_img,bf_blobs);

    /* Display result */
    flip_image(imCapt,imDisp);
    paste_image(imDisp,bf_blobs,width);
    
    if(kb_lowresfl) {
      upsample_image(imDisp,imDisp2);
      cvShowImage(window1, imDisp2);
    } else {
      cvShowImage(window1, imDisp);
    }

    ikey=cvWaitKey(5); /* Needed for highgui event processing */
    if(ikey>0) {
      key=(char)ikey;

      if(key == 'n') {
	kb_n8flag=1-kb_n8flag;
	printf("n8flag=%d\n",kb_n8flag);
      }
      
      if(key == 'd') {
	kb_normfl=1-kb_normfl;
	printf("normfl=%d\n",kb_normfl);
      }
      
      if(key == 'b') {
	kb_blurfl=1-kb_blurfl;
	printf("blurfl=%d\n",kb_blurfl);
      }
      if(key =='c') {
	kb_cspace=1-kb_cspace; /* Toggle colourspace */
	printf("cspace=%d\n",kb_cspace);
      }
      
      if(key == 'r') {
	tb_margin=DEF_MARGIN; /* Reset to default */
	tb1_val=tb_margin/TB1_SCALE*TB1_MAX;
	cvSetTrackbarPos(trackbar1,window1,tb1_val);
	tb2_val=tb_timesteps=DEF_TIMESTEPS;  /* Reset to default */
	cvSetTrackbarPos(trackbar2,window1,tb2_val);
	kb_n8flag=DEF_N8FLAG;
	kb_normfl=DEF_NORMFL;
	kb_blurfl=DEF_BLURFL;
	printf("timesteps =%d\n",tb_timesteps);
	printf("min_margin=%g\n",tb_margin);
	printf("n8flag=%d\n",kb_n8flag);
	printf("normfl=%d\n",kb_normfl);
	printf("blurfl=%d\n",kb_blurfl);
	cvWaitKey(1); 
	/* printf("tb1:%d\n",cvGetTrackbarPos(trackbar1,window1)); */
	/* printf("tb2:%d\n",cvGetTrackbarPos(trackbar2,window1)); */
      }
      if(key == 'h') {
	kb_lowresfl=1-kb_lowresfl; /* Toggle resolution */
	if(kb_lowresfl) {
	  width=320;
	  height=240;
	} else {
	  width=640;
	  height=480;
	}
	cvReleaseCapture(&cam);
	cam = cvCreateCameraCapture(0);
	/* Set size of image */
	cvSetCaptureProperty(cam, CV_CAP_PROP_FRAME_WIDTH, width);
	cvSetCaptureProperty(cam, CV_CAP_PROP_FRAME_HEIGHT, height);
	/* Free image buffers */
	cvReleaseImage( &imDisp );
	buffer_free(bf_img);
	bbuffer_free(bf_blobs);
	/* Allocate image buffers */
	imDisp =  cvCreateImage( cvSize(width*2,height), IPL_DEPTH_8U, 3);
	bf_img = buffer_new(height,width,3);
	bf_blobs = bbuffer_new(height,width,3);
      }
    }

    timewin_addtime(tw);
    fno++;
    if((fno%DEF_TIMESIZE)==0) printf("Frame rate %g Hz\n",timewin_rate(tw));
  }
  /* Clean up */
  timewin_free(tw);
  buffer_free(bf_img);
  bbuffer_free(bf_blobs);
  cvDestroyAllWindows();
  cvReleaseCapture(&cam);
  return 0;
}
