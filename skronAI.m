function M=skronAI(A)
%
%  Compute the symetric kronnecker product of A and identity
%
s=size(A);
ir=s(1);
ic=s(2);
if ir~=ic
	error('Matrice is not square');
else;
	Q=vecPsvec(ic);
	C=skronaipart(A,ic);
	M=Q'*C*Q;
end

