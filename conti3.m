
trial=3;
bigrep=1;
R=3;
V=50;
n=100
pp=0.3;
sigma_error=1;


figamma = fopen(strcat('cerror',num2str(trial),'.txt'),'w');
for rep=1:bigrep
fic = fopen(strcat('raw',num2str(bigrep),'.txt'),'w');
fia = fopen(strcat('app',num2str(bigrep),'.txt'),'w');

    E=zeros(2,4);

rand('seed',rep+bigrep*(trial-1));
randn('seed',rep+bigrep*(trial-1));

    B = zeros(V,R);
    for i=1:R
        B(:,i)=randn(V,1);
    end

    Lambda=zeros(R,R,n);
    L=zeros(V,V,n);
    tL=zeros(V,V,n);
x=randn(n,1)+0.5;
%x=binornd(1,0.5,[n,1]);
for i=1:n
        A=zeros(R);

        A(triu(true(size(A)))) = randn(R*(R+1)/2,1)*1/sqrt(2)+1;
A(1,2)=A(1,2)+x(i,1)*4;
A(2,3)=A(2,3)+x(i,1)*4;
        Lamb=A+A.'-2*diag(diag(A));

        Lamb(1:(R+1):end)=(normrnd(1,1,R,1));
        Lambda(:,:,i)=Lamb;

        a=normrnd(1,sigma_error/sqrt(2),V*(V+1)/2,1);
        AA=zeros(V);
        AA(triu(true(size(AA)))) = a;
        noise=AA+AA.'-2*diag(diag(AA));
        noise(1:(V+1):end)=normrnd(0,sigma_error,V,1);
      %  Z=diag([ones(1,pp*R) zeros(1,(1-pp)*R)]);
        tL(:,:,i)=B*Lambda(:,:,i)*B';
        L(:,:,i)=B*Lambda(:,:,i)*B'+noise;
    end

%fprintf(fic,[repmat('%4.4f\t',1,2),'\n'],L(:,:,1));
%[nB,nLambda,nL,nGamma1,nGamma2,dic,bic,aic,error] =HGSC(L,([ones(n,1),x ]),3,1,500,100,rep+bigrep*(trial-1),B,Lambda);

MCMCpara.Niter=30;
MCMCpara.burnin=0;
MCMCpara.B=zeros(V,R);
MCMCpara.Lambda=zeros(R,R,V);
MCMCpara.Gamma=zeros(2,R*(R+1)/2);

MCMCpara.sigma=1;
MCMCpara.sig_gam=1;
MCMCpara.b1=10^(-2);
MCMCpara.b2=10^(-2);
MCMCpara.c1=10^(-2);
MCMCpara.c2=10^(-2);
MCMCpara.va=1/2;
MCMCpara.vb=1/2;
slice.miter=10
slice.burnin=0;
slice.width=10;
Output= PX_BLGRM(L,[ones(100,1),x],3,MCMCpara,slice);

%fprintf(fia,[repmat('%4.4f\t',1,2),'\n'],nL(:,:,1));
disp('error');
nGamma1=zeros(3);
nGamma1(tril(true(size(nGamma1)))) = Output.nGamma(1,:);
nGamma1=nGamma1+nGamma1.'-diag(diag(nGamma1));

nGamma2=zeros(3);
nGamma2(tril(true(size(nGamma2)))) = Output.nGamma(2,:);
nGamma2=nGamma2+nGamma2.'-diag(diag(nGamma2));

ner1=norm(B*[0,4,0;4,0,4;0,4,0]*B'-Output.nB*nGamma2*Output.nB.','fro')/norm(B*[0,4,0;4,0,4;0,4,0]*B','fro');
ner2=norm(B*ones(3,3)*B'-Output.nB*nGamma1*Output.nB.','fro')/norm(B*ones(3,3)*B','fro');
disp(ner1);
E(1,rep)=ner2;
E(2,rep)=ner1;


%for rr=1:(pp*R*2)

%    [nB,nLambda,nL,dic,bic,aic]=GSC0(L,rr,5000,2000);

 %   a=0;
 %   for k=1:n
 %       a=a+(norm(tL(:,:,k)-nL(:,:,k),'fro')/norm(tL(:,:,k),'fro'));
 %   end
 %   disp(a/n);E(1,1)=a/n;
  %  E(1,2)=dic; E(1,3)=bic; E(1,4)=aic;
%fprintf(figamma,[repmat('%4.4f\t',1,4),'\n'],E);
%end
fprintf(figamma,[repmat('%4.4f\t',1,2),'\n'], [ner1;ner2]);
end

