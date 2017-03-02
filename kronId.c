/* Compute the Kronecker product of (double) matrix A and identity 
 *
 * Input : 
 *
 * pointer to the full matrix A
 * number of row of Id matrix ir
 * number of column of Id matrix ic
 * number of row of A matrix sar
 * number of column of A matrix sac
 *
 * Output : 
 *
 * pointer to sparse matrix format int* I (row), int* J (column), double K
 * (value)
 *
 * */
void kronAI(double *A, unsigned int ir, unsigned int ic, unsigned int sar, unsigned int sac, int *I, int *J, double *K)
{
    int k,p,pos1,pos2,mi,n;

	// find min dimension of I
    mi = (ir<ic) ? ir :  ic;

    // compute A kron I
    for (k=0;k<sar;k++)
    {
		pos1=k*ir;
		for (p=0;p<sac;p++)
	  	{
	    	pos2=p*ic;
	    	for (n=0;n<mi;n++)
	      	{
				*I=pos1+n;
				*J=pos2+n;
				*K=A[k+p*sar];
				I++;
				J++;
				K++;
	      	}
	  	}
    }
}



/* Compute the kronecker product of identity and (double) matrix A 
 *
 * pointer to the full matrix A
 * number of row of Id matrix ir
 * number of column of Id matrix ic
 * number of row of A matrix sar
 * number of column of A matrix sac
 *
 * Output : 
 *
 * pointer to sparse matrix format int* I (row), int* J (column), double K
 * (value)
 *
 * */
void kronIA(double *A, unsigned int ir, unsigned int ic, unsigned int sar, unsigned int sac, int *I, int *J, double *K)
{
    int k,p,pos1,pos2,mi,n;

    // find min dimension of I
    mi = (ir<ic) ? ir :  ic;

    // compute I kron A
    for (k=0;k<mi;k++)
    {
		pos1=k*sar;
		pos2=k*sac;
		for (p=0;p<sar;p++)
	  	{
	    	for (n=0;n<sac;n++)
	      	{
				*I=pos1+p;
				*J=pos2+n;
				*K=A[p+n*sar];
				I++;
				J++;
				K++;
	      	}
	  	}
    }
}

/* Compute (double) A skron I, symmetric kronecker procduct 
 *
 * pointer to the full matrix A
 * number of row of Id matrix ic (symmetric)
 *
 * Output : 
 *
 * pointer to full (double) matrix C
 * */
void skronaipart(double *A, unsigned int ic, double *C)
{
    int k,p,pos1,pos2,ncc,n;

    // Compute size of output matrix C
    ncc=ic*ic;
    // compute A skron I -> to use with vecPsvec
    for (k=0;k<ic;k++)
    {
		pos1=k*ic;
		for (p=0;p<ic;p++)
	  	{
	    	pos2=p*ic;
	    	for (n=0;n<ic;n++)
	      	{
				C[pos1*ncc+n*ncc+pos2+n]=A[k*ic+p]/2;
	      	}
	  	}
    }
   	for (k=0;k<ic;k++)
    {
		pos1=k*ic;
		pos2=k*ic;
		for (p=0;p<ic;p++)
	  	{
	    	for (n=0;n<ic;n++)
	      	{
				C[pos1*ncc+p*ncc+pos2+n]+=A[p*ic+n]/2;
	      	}
	  	}
    }
}

/* svec : symmetric vectorisation 
 *
 * Inputs:
 *
 * int n : the size of the symmetric matrix
 * double *M : pointer to the matrix
 *
 * Outputs:
 *
 * double *v : pointer to the output vector v
 *
 * */
void svecmex(int n,double *M, double *v)
{
    int i,j,pos;
    pos=0;
    for (i=0;i<n;i++)
    {
        for (j=i;j<n;j++)
        {
            if (i==j)
            {
                v[pos]=M[n*i+j];
            }
            else
            {
                v[pos]=M[n*i+j]*1.4142135623730950488;
            }
            pos++; 
        }
    }
}

/* sveci : inverse of symmetric vectorisation 
 *
 * Inputs:
 *
 * int n : the size of the symmetric matrix
 * double *v : pointer to the output vector v
 *
 * Outputs:
 *
 * double *M : pointer to the matrix
 *
 * */
void svecimex(int n,double *v, double *M)
{
    int i,j,pos;
    pos=0;
    for (i=0;i<n;i++)
    {
        for (j=i;j<n;j++)
        {
            if (i==j)
            {
                M[n*i+j]=v[pos];
            }
            else
            {
                M[n*i+j]=v[pos]/1.4142135623730950488;
                M[n*j+i]=v[pos]/1.4142135623730950488;                
            }
            pos++; 
        }
    }
}


/* vecPsvec
 *
 * Compute the P matrix s.t. 
 *
 *   vec(X)= P * svec(X)
 *
 * where X is a symmetric matrix of size n
 *
 * INPUT:
 * 
 * - (int) n
 *
 * OUTPUT: 
 * 
 * - pointers to sparse matrix format of P.
 *   P is of size n*n x m, and memory of sparse matrix is n*n
 *
 */
 void vecPsvec(unsigned int n, int *I, int *J, double *K)
{
	int i,seq,test,co,co2,co3,co4=1;
	short sw=0;
	seq=0;
	test=0;
	co=n;
	for (i=0;i<n*n;i++)
	{
		if (i==test)
		{
			// diagonal 

			// reset seq after lower diagonal part	
			if (i==0)
			{
				co3=0;
			}
			else
			{
				co3+=co+1;
			}
			seq=co3;

			I[i]=i;
	    	J[i]=seq;
			K[i]=1.0;
			test+=n+1;
			seq++;
			co--;
			sw=0;
		}
		else
		{
			if (i==test-n+co) 
			{
				sw=1;
				seq=co4;
				co4++;
				co2=n-1;
			}
			// diagonal 
			if (sw==0)
			{
				// upper diagonal
				I[i]=i;
	    		J[i]=seq;
				K[i]=1/1.4142135623730950488;
				seq++;
			}
			else
			{
				// lower diagonal
				I[i]=i;
	    		J[i]=seq;
				K[i]=1/1.4142135623730950488;
				seq+=co2;
				co2--;
			}
		}
	}
}
