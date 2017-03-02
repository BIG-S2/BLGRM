function S=sk(b)
k=size(b,2);
R=k*(k+1)/2;
S=eye(R);
l=0;
for i=1:k
for j=i:k
    l=l+1;
S(l,l)=b(i)*b(j);
end
end