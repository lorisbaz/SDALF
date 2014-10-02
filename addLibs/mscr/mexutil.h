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

#ifndef _MEXUTIL
#define _MEXUTIL

#include "image_buffer.h"

buffer *buffer_encapsulate(const mxArray *arg);
ibuffer *ibuffer_encapsulate(const mxArray *arg);

#define SSTR_SIZE 30
int type_check(char *pname,const mxArray *arg,mxClassID classid,int mdims,
	       int *msize);

double get_field_value(const mxArray *matlab_struct,char *field_name,
		       double defval);

#endif /* _MEXUTIL */
