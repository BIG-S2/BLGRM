function S=vecPsvec(n)
[I,J,K]=vecPsvecmex(n);
S=sparse(I,J,K);
