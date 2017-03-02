#include "mex.h"

/* prototype */
void svecmex(int n,double *M, double *v);

/* Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    double *A,*Ip;
    int ir,nzmax;
    /* Check Input-Output */
    if (nrhs!=1)
    {
        mexErrMsgTxt("One input required: size n, and matrix A (real, symmetric) ");
    }
    else if (nlhs>1)
    {
        mexErrMsgTxt("Too many outputs arguments.");
    }
	/* read input data */
    A=mxGetPr(prhs[0]);
    ir=(int)mxGetM(prhs[0]);
   
	nzmax = ((ir+1)*ir)/2;
    plhs[0]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
	
    Ip=mxGetPr(plhs[0]);

	svecmex(ir,A, Ip);
  }
