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
#include <string.h>
#include "mex.h"
#include "image_buffer.h"
#include "mexutil.h"

/*
 *
 */
buffer *buffer_encapsulate(const mxArray *arg) {
const int *dims;
buffer *bf_ret;
   dims = mxGetDimensions(arg); 
   switch(mxGetNumberOfDimensions(arg)) {
   case 0:
     bf_ret = buffer_new0((fpnum *)mxGetPr(arg),1,1,1);
     break;
   case 1:
     bf_ret = buffer_new0((fpnum *)mxGetPr(arg),dims[0],1,1);
     break; 
   case 2:
     bf_ret = buffer_new0((fpnum *)mxGetPr(arg),dims[0],dims[1],1);
     break;
   default:
     bf_ret = buffer_new0((fpnum *)mxGetPr(arg),dims[0],dims[1],dims[2]);
   }
   return bf_ret;
}
/*
 *
 */
ibuffer *ibuffer_encapsulate(const mxArray *arg) {
const int *dims;
ibuffer *bf_ret;
   dims = mxGetDimensions(arg); 
   switch(mxGetNumberOfDimensions(arg)) {
   case 0:
     bf_ret = ibuffer_new0((int *)mxGetPr(arg),1,1,1);
     break;
   case 1:
     bf_ret = ibuffer_new0((int *)mxGetPr(arg),dims[0],1,1);
     break; 
   case 2:
     bf_ret = ibuffer_new0((int *)mxGetPr(arg),dims[0],dims[1],1);
     break;
   default:
     bf_ret = ibuffer_new0((int *)mxGetPr(arg),dims[0],dims[1],dims[2]);
   }
   return bf_ret;
}
/*
** Function to check wether an argument is
** of proper class and size
** -1 as size means any size
*/
#define SSTR_SIZE 30
int type_check(char *pname,const mxArray *arg,mxClassID classid,int mdims,int *msize)
{
const int *dims;
 const char *classidl[]={"double","int","uint8","other"};
 char sstr[SSTR_SIZE+1];
 int k,retval;
  dims = mxGetDimensions(arg);
  retval=1;
  if(mxGetClassID(arg) != classid) {
    switch(classid) {
    case mxDOUBLE_CLASS:
      k=0;break;
    case mxINT32_CLASS:
      k=1;break;
    case mxUINT8_CLASS:
      k=2;break;
    default:
      k=3;
    }
   mexPrintf("Type mismatch: <%s> should be of class %s.\n",pname,classidl[k]);
    retval=0;
  }
  if(mxGetNumberOfDimensions(arg) != mdims) {
    mexPrintf("Type mismatch: <%s> should be %dD.\n",pname,mdims);
     retval=0;
  }
  if(retval) {
    for(k=0;k<mdims;k++) {
      if((msize[k] != -1) && (dims[k] != msize[k])) retval=0;
    }
    if(!retval) {
      strncpy(sstr,"proper size",SSTR_SIZE);
      if((mdims==2)&&(msize[0]==1)&&(msize[1]==1))
	strncpy(sstr,"scalar",SSTR_SIZE);
      if((mdims==2)&&(msize[0]>0)&&(msize[1]>0))
	sprintf(sstr,"%dx%d",msize[0],msize[1]);
      if((mdims==3)&&(msize[0]>0)&&(msize[1]>0)&&(msize[2]>0))
	sprintf(sstr,"%dx%dx%d",msize[0],msize[1],msize[2]);
      mexPrintf("Type mismatch: <%s> should be %s.\n",pname,sstr);
    }
  }
  return retval;
}
/* 
 * Get scalar value from a field in a matlab struct
 */
double get_field_value(const mxArray *matlab_struct,char *field_name,double def_val) {
double value;
int field_no;
  field_no = mxGetFieldNumber(matlab_struct,field_name);
  if(field_no>=0) {
    value = mxGetScalar(mxGetFieldByNumber(matlab_struct,0,field_no));
  } else {
    value = def_val;
  }
  return value;
}
