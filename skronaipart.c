#include "mex.h"

/* prototype */
void skronaipart(double *A, unsigned int ic, double *C);

/* Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    double *A,*Ip;
    int ic,nzmax;
    /* Check Input-Output */
    if (nrhs!=2)
    {
        mexErrMsgTxt("two inputs required:  matrix A (real, symmetric) ");
    }
    else if (nlhs>1)
    {
        mexErrMsgTxt("Too many outputs arguments.");
    }
	/* read input data */
    A=mxGetPr(prhs[0]);
    ic=(int)mxGetScalar(prhs[1]);
  
	/* memory allocation */
	nzmax = ic*ic;	
    plhs[0]=mxCreateDoubleMatrix(nzmax,nzmax,mxREAL);
	
    Ip=mxGetPr(plhs[0]);

	skronaipart(A, ic, Ip);
  }
