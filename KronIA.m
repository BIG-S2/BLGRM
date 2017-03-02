function C=KronIA(A,ir,ic)
% 
% function C=KronIA(A,ir,ic)
%
% Compute the kronecker product I \otime A, of identity matrix I
% and matrix A. ir is the number of row of I, ic is the number of
% column of C. C is a sparse matrix.
%
[I,J,S]=kronIAmex(A,ir,ic);
C=sparse(I,J,S,ir*size(A,1),ic*size(A,2));

