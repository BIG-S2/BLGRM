#include "mex.h"
#include <math.h>

/* prototype */
void svecimex(int n,double *v, double *M);

/* Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    double *A,*Ip;
    int ir,nzmax;
    /* Check Input-Output */
    if (nrhs!=1)
    {
        mexErrMsgTxt("One input required: matrix A (real, symmetric) ");
    }
    else if (nlhs>1)
    {
        mexErrMsgTxt("Too many outputs arguments.");
    }
	/* read input data */
    A=mxGetPr(prhs[0]);
    ir=(int)mxGetM(prhs[0]);
  
    nzmax=(sqrt(8*ir+1)-1)/2;
    plhs[0]=mxCreateDoubleMatrix(nzmax,nzmax,mxREAL);
	
    Ip=mxGetPr(plhs[0]);

	svecimex(nzmax,A, Ip);
  }
