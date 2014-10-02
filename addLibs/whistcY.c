/* file:        whistcY.c
** description: MEX weighted histc function.
**/


/** @file
 ** @brief HISTC MEX function implementation.
 **/

#include"mex.h"
#include<stdlib.h>
#include<math.h>

/** WHISTC(X,W,EDGES) 
 **/
#define min(a,b) ((a<b)?a:b)
#define max(a,b) ((a<b)?b:a)


/** @brief MEX driver.
 ** @param nout number of MATLAB output arguments.
 ** @param out MATLAB output arguments.
 ** @param nin number of MATLAB input arguments.
 ** @param in MATLAB input arguments.
 **/  
void
mexFunction(int nout, mxArray *out[], 
            int nin, const mxArray *in[])
{
  int M, N, NE ;
  double* Xpt ;
  double* Wpt ; 
  double* EDGESpt ;
  double* RESpt ;
  enum {X=0, W, EDGES} ;

  /** -----------------------------------------------------------------
   **                                               Check the arguments
   ** -------------------------------------------------------------- */
//   if (nin != 3) {
//     mexErrMsgTxt("Three arguments required.");
//   } else if (nout > 1) {
//     mexErrMsgTxt("Too many output arguments.");
//   }

//   if (!mxIsDouble(in[X]) || 
//       !mxIsDouble(in[W]) ||
//       !mxIsDouble(in[W])) {
//     mexErrMsgTxt("The arguments must be real matrices.") ;
//   }

  M = mxGetM(in[X]) ;
  N = mxGetN(in[X]) ;
//   if( M != mxGetM(in[W]) ||
//       N != mxGetN(in[W]) ) {
//     mexErrMsgTxt("X and W must have the same dimension.") ;
//   }
// 
//   if(min( mxGetM(in[EDGES]), mxGetN(in[EDGES]) ) != 1) {
//     mexErrMsgTxt("EDGES must be a vector.") ;
//   }
  
  NE = max(mxGetM(in[EDGES]), mxGetN(in[EDGES])) ;

//   if(NE < 2) {
//     mexErrMsgTxt("At least two edges are required.\n") ;
//   }
  
  Xpt = mxGetPr(in[X]) ;
  Wpt = mxGetPr(in[W]) ;
  EDGESpt = mxGetPr(in[EDGES]) ;

//   {
//     double x = EDGESpt[0] ;
//     int j ;
//     int ok = 1 ; 
//     for(j = 1 ; j < NE ; ++j) {
//       ok &= (x < EDGESpt[j]) ;
//       x = EDGESpt[j] ;
//     }
//     if(!ok) mexErrMsgTxt("EDGES must be increasing.") ;
//   }

  /*if(nout == 0) return ;*/

  /* If the input is a vector, make it a column. */
  if(M == 1) {
    M = N ; 
    N = 1 ;
  }

  /* Alloc the result. */
  out[0] = mxCreateDoubleMatrix(NE, 1, mxREAL) ;
  RESpt = mxGetPr(out[0]) ; 

  /** -----------------------------------------------------------------
   **                                                        Do the job
   ** -------------------------------------------------------------- */
	{  
     int i = 0;
	  for(; i < M; i++){
			int c = 0;
//             mexPrintf("%f    ",Xpt[i]);
			while(c < NE){
				if((c==NE-1) || (Xpt[i] < EDGESpt[c+1] && Xpt[i] >= EDGESpt[c])){
					RESpt[c] += Wpt[i];
					c = NE+1; // exit from while
				}
				else
					c++;
			  }
	  }
     
  }
  return ;
}


