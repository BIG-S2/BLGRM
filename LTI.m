 function A=LTI(a,A1,i,V,R)
    A=zeros(V,R);

A(i:V,i)=a.';
A1(:,i)=0;
A=A+A1;
    end