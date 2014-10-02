#
# Makefile for MSCR feature detector
#
VERSION = 1.01

# Add preprocessor flags here
DEFS = -DDOUBLE_FLOATS -DINTEGER_EDGES # For Matlab targets
#DEFS = -DINTEGER_EDGES  # Use this instead for faster opencv target

# GCC compiler
CC= gcc

# Location of OpenCV
#OPENCV = /usr/local/include/opencv/
#OPENCV = /cs/public/generic/lib/pkg/opencv-1.0.0/include/opencv/
OPENCV = /usr/include/opencv/

# Location of Matlab
MATLAB=/cs/local/generic/lib/pkg/matlab-7.0.1
#MATLAB=/Applications/MATLAB704

# Linker paths
#LDPATHS= -L/usr/local/lib/
#LDPATHS = -L/cs/public/lib/
#LDPATHS = -L/cs/public/lib/pkg/opencv-1.0.0/lib/
LDPATHS = -L/usr/lib/

# Linker flags for opencv targets
LDFLAGS1=$(LDPATHS) -lm -lc -lhighgui -lcv -lcxcore

# Linker flags for matlab targets
LDFLAGS2=$(LDPATHS) -lm -lc

# Matlab mex compiler script
MEX=${MATLAB}/bin/mex
MEXFLAGS=-I${MATLAB}/extern/include/ ${DEFS} -v

# Include paths
INCPATHS =  -I${OPENCV} -I${MATLAB}/extern/include/

#CFLAGS= -O3 -Wall -ffast-math -fforce-addr ${DEFS} -mtune=prescott ${INCPATHS}
CFLAGS= -g -Wall  ${DEFS}  ${INCPATHS}

# Object files needed for building mscr_opencv app
OPENCVOBJ = mscr_opencv.o measure_time.o image_buffer.o msr_util.o visualise.o

# Object files needed for detect_mscrmex app
MSCRMEXOBJ = image_buffer.o msr_util.o mexutil.o

# Object files needed for draw_blobs app
DRAWMEXOBJ = image_buffer.o visualise.o mexutil.o

# C source
CFILES = mscr_opencv.c detect_mscr.c measure_time.c msr_util.c image_buffer.c visualise.c mexutil.c draw_blobs.c

# H files
HFILES = chi_table.h image_buffer.h measure_time.h msr_util.h visualise.h mexutil.h

# M files
MFILES = matlab_demo.m detect_mscr.m draw_blobs.m blobs_inside_image.m

# Auxiliary files
AUXFILES = Makefile README COPYING ChangeLog.txt puzzle.png

# Full distribution
DISTR_FILES = ${CFILES} ${HFILES} ${MFILES}  ${AUXFILES}

DISTR_NAME = mscr_${VERSION}

TARGETS = detect_mscr draw_blobs mscr_opencv
all:	${TARGETS}

clean:
	rm *.o *.mexglx ${TARGETS}

distr:
	tar zcvf ${DISTR_NAME}.tar.gz $(DISTR_FILES)
	zip ${DISTR_NAME}.zip $(DISTR_FILES)

mscr_opencv:	$(OPENCVOBJ)
	$(CC) $(CFLAGS) -o $@ $(OPENCVOBJ) $(LDFLAGS1)

detect_mscr:	$(MSCRMEXOBJ) detect_mscr.c
	$(MEX) $(MEXFLAGS) detect_mscr.c -output $@ $(MSCRMEXOBJ) $(LDFLAGS2)

draw_blobs:	$(DRAWMEXOBJ) draw_blobs.c
	$(MEX) $(MEXFLAGS) draw_blobs.c -output $@ $(DRAWMEXOBJ) $(LDFLAGS2)

%.o:	%.c
	$(CC) $(CFLAGS) -c $<
