function C=KronAI(A,ir,ic)
% 
% function C=KronAI(A,ir,ic)
%
% Compute the kronecker product A \otime I, of identity matrix I
% and matrix A. ir is the number of row of I, ic is the number of
% column of C. C is a sparse matrix.
%
[I,J,S]=kronAImex(A,ir,ic);
C=sparse(I,J,S,ir*size(A,1),ic*size(A,2));

