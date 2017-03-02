#include "mex.h"
#include <stdlib.h>

/* prototype */
void kronAI(double *A, unsigned int ir, unsigned int ic, unsigned int sar, unsigned int sac, int *I, int *J, double *K);

/* Gateway routine */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,const mxArray *prhs[])
{
    double *A,*K,*Ip,*Jp,*Kp;
	int	*I, *J;
    int ir,ic,sa1,sa2,mi,nzmax,nrc,ncc;
    /* Check Input-Output */
    if (nrhs!=3)
    {
        mexErrMsgTxt("three inputs required: A, i_row and i_column");
    }
    else if (nlhs>3)
    {
        mexErrMsgTxt("Too many outputs arguments.");
    }
	/* read input data */
    A=mxGetPr(prhs[0]);
    ir=(int)mxGetScalar(prhs[1]);
    ic=(int)mxGetScalar(prhs[2]);
    sa1=mxGetM(prhs[0]);
    sa2=mxGetN(prhs[0]);
   
   	/* Create output sparse matrix */
	
	/* compute size of the output matrix  */
	nrc=ir*sa1;
    ncc=ic*sa2;

	/* find min dimension of I */
    mi = (ir<ic) ? ir :  ic;
	
	/* allocation size  */
	nzmax = mi*sa1*sa2;
    plhs[0]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[1]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
    plhs[2]=mxCreateDoubleMatrix(nzmax,1,mxREAL);
	
    Ip=mxGetPr(plhs[0]);
    Jp=mxGetPr(plhs[1]);
    Kp=mxGetPr(plhs[2]);
	
	/* pointer allocation */
	I = malloc(nzmax*sizeof(int));
	J = malloc(nzmax*sizeof(int));
	K = malloc(nzmax*sizeof(double));

	kronAI(A, ir, ic, sa1, sa2, I, J, K);

	for (ir=0;ir<nzmax;ir++)
	{
		Ip[ir] = (double) I[ir]+1;
		Jp[ir] = (double) J[ir]+1;
		Kp[ir] = K[ir];
	}

	/* free memory */
	free(I);
	free(J);
	free(K);
  }
