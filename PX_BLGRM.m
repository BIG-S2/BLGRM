function Output = PX_BLGRM (L,X,R,MCMCpara,slice)
%----------------------------------------------------------------------------------
%		PX-Bayesian Low-rank Graph Regression Model	   							  %
%----------------------------------------------------------------------------------
%
% Decomposition model: L_i = B*Lambda_i*B + E		   
% Regression model: Lambda_i = X*Gamma + e 		   
%----------------------------------------------------------------------------------
%
% Data input:
%	L: V x V x n, input data array. 
%	   A three-dimensional array with n symmetric V x V matrices.
%	X: n x p, covariate matrix for the regression model.
%	n: the number of subjects.
%	V: the number of edges or the dimension of input matrix.	
%	R: the prespecified number of eigenvectors. We use BIC to determine R.	
%----------------------------------------------------------------------------------
%
% MCMC setting input:
%	MCMCpara: input structure of the total number of MCMC iterations, burn-in, 
%			  initial values, hyperparameters for MCMC.
%	 	
%		MCMCpara.Niter: the total number of MCMC iterations.
%		MCMCpara.burnin: the number of burn-in
%		MCMCpara.B: the initial value for B matrix.
%		MCMCpara.Lambda: the initial value for Lambda matrix.
%		MCMCpara.Gamma: the initial value for Gamma matrix.
%		MCMCpara.sigma: the initial value for sigma (decomposition model variance).
%		MCMCpara.sig_gam: the initial value for sig_gam (regression model variance).
%		MCMCpara.b1: the initial value for hyperparameter for priors of sigma.
%		MCMCpara.b2: the initial value for hyperparameter for priors of sigma.
%		MCMCpara.c1: the initial value for hyperparameter for priors of sig_gam.
%		MCMCpara.c2: the initial value for hyperparameter for priors of sig_gam.
%		MCMCpara.va: the initial value for parameter for priors of Psi.
%		MCMCpara.vb: the initial value for parameter for priors of Psi.
%----------------------------------------------------------------------------------
%
% Slice sampling setting input:
%	slice: input structure of tuning parameters for slice sampling for B.
%		slice.miter: the maximum number of iterations for tuning. 
%		slice.burnin: the number of burnin during the tuing iterations. 
%		slice.width: the width of slice sampling
%----------------------------------------------------------------------------------
	
tic;
devhat=0;
sigmahat=0;

V=size(L,1);
n=size(L,3);
sigma0=1;	% We fix the variance of Lambda prior for identifiability
rndd=ones(R,1)*2;
Psi=diag(rndd);	% initialization of Psi
lambb=gamrnd(1,1/2,[R,1]);

%%% Initialization 
Niter=MCMCpara.Niter;
burnin=MCMCpara.burnin;
B=MCMCpara.B;
Lambda=MCMCpara.Lambda;
Gamma=MCMCpara.Gamma;

sigma=MCMCpara.sigma;
sig_gam=eye(R*(R+1)/2,R*(R+1)/2)*MCMCpara.sig_gam;
b1=MCMCpara.b1;
b2=MCMCpara.b2;
c1=MCMCpara.c1;
c2=MCMCpara.c2;
va=MCMCpara.va;
vb=MCMCpara.vb;

phi=ones(V,R);
for i=1:(R-1)
	for j=(i+1):R
		phi(i,j)=0;;
	end
end

tau=ones(R,1);

%%% Memory allocation 
aLambda=zeros(R,R,n);
sLambda=zeros(R,R,n);
sB=zeros(V,R);
sL=zeros(V,V,n);
nL=zeros(V,V,n);
LL=zeros(V,V,n);
Delta=zeros(R,R,n);
p=size(X,2);
sGamma=zeros(p,R*(R+1)/2);
mug=zeros(p*R*(R+1)/2,1);          
sg=zeros(p*R*(R+1)/2,p*R*(R+1)/2);

gamtau=ones(p*R*(R+1)/2,1);
rnd=ones(R,1);
lamb=randn(R*(R+1)/2,1)/sqrt(sigma0)/sqrt(2);
%%% D matrix calculation (See eq.(6)) 
a=vecPsvec(R);
gg=sum(a,1);
D=diag(gg);


for iter=1:Niter

    tic;
    U=B.'*B;
	cov=D*skron(U,U)*D*sigma+sigma0*D*skron(Psi,Psi)*D;
	cov=tril(cov)+tril(cov).'-diag(diag(cov));	% Force to be symmetric
    Cov=inv(cov);
    Cov=tril(Cov)+tril(Cov).'-diag(diag(Cov));
    C=cov\D;

%%% Update Lambda    
    for j=1:n
        W=B.'*L(:,:,j)*B*sigma+Psi*Delta(:,:,j)*Psi*sigma0;
		M=C*svecmex(W);
		rr=0;
		for i=1:R
			gg=lamb;
			gg(rr+(1:i))=[];
			ss=Cov(rr+(1:i),:);
			ss(:,rr+(1:i))=[];
			lamb(rr+(1:i))=mvnrnd(M(rr+(1:i))-ss*gg, Cov(rr+(1:i),rr+(1:i))).';
			rr=rr+i;
		end
        A=zeros(R);
        A(tril(true(size(A))))=lamb;
      	Lambda(:,:,j)=A+A.'-diag(diag(A));
		%% mug and sg will be used for update of Gamma
		mug=mug+reshape(X(j,:)'*svecmex(Psi*Lambda(:,:,j)*Psi)'*(D),[p*R*(R+1)/2,1]);	
        sg=sg+kron(D*skron(Psi,Psi)*D,X(j,:)'*X(j,:));
    end
	
%%% Update Gamma    
	S_gam=inv(sigma0*sg+kron(sig_gam,eye(p,p)));
  	S_gam=tril(S_gam)+tril(S_gam).'-diag(diag(S_gam));	% Force to be symmetric
    mu_gam=sigma0*S_gam*mug;
    for i=1:(R*(R+1)/2)
		rr=(p*(i-1)+1):(p*(i-1)+p);
		gg=reshape(Gamma,[p*R*(R+1)/2,1]);
        gg(rr)=[];
        ss=S_gam(rr,:);
        ss(:,rr)=[];
        Gamma(:,i)=mvnrnd(mu_gam(rr)-ss*gg, S_gam(rr,rr)).';
	end
                                  
%%% Calculate Delta    
    for j=1:n
		AD=zeros(R);
        AD(tril(true(size(AD))))=X(j,:)*Gamma;
        Delta(:,:,j)=AD+AD.'-diag(diag(AD));
    end

	
%%% Update B using slice sampling    
    lams=zeros(R^2,R^2);
    laml=zeros(V*R,V*R);
    
    for  j=1:n
        lams=lams+kron(Lambda(:,:,j),Lambda(:,:,j));
        laml=laml+kron(Lambda(:,:,j),L(:,:,j));
    end
	
	JJ=phi;
	VV0=diag(JJ(tril(true(size(JJ)))));
	why0=1;
	for i=1:R
		why=i*V-i*(i-1)/2;
    	logpdf = @(b1) 1/2*n*V*(V+1)*log(sigma)-sigma/2*(-2*reshape(LTI(b1,B,i,V,R),[V*R,1]).'*laml*reshape(LTI(b1,B,i,V,R),[V*R,1])+reshape(LTI(b1,B,i,V,R).'*LTI(b1,B,i,V,R),[R^2,1]).'*lams*reshape(LTI(b1,B,i,V,R).'*LTI(b1,B,i,V,R),[R^2,1]))-1/2*b1*VV0(why0:why,why0:why)*b1.';
		
		if iter<=slice.miter
			rnd=slicesample(B(i:V,i).',1,'logpdf',logpdf,'burnin', slice.burnin,'width',slice.width);
			BB00=B;
		else
			rnd=slicesample(B(i:V,i).',1,'logpdf',logpdf,'burnin', 0,'width',slice.width);
		end
			
		B(i:V,i)=rnd;
		why0=why+1;
	end
    
  
	    
%%%Update hyperparameters
	for i=1:V
		for j=1:min(i,R)
			phi(i,j)=random('inversegaussian',sqrt(lambb(j)/B(i,j)^2),lambb(j));
		end
	end
	
	for i=1:R
    	bv=B(i:V,i);
    	tau(i)=gamrnd((V-i+1)/2,1/(bv.'*diag(phi(i:V,i))*bv/2));
	end
	
	for i=1:R
    	lambb(i)=gamrnd(1/2+(V-i+1),1/(sum(((1./phi(i:V,i))))/2+1));
	end
	
	BMM=zeros(R*(R+1)/2,R*(R+1)/2);
	for j=1:n
    	sl=svecmex(Lambda(:,:,j)-Delta(:,:,j));
    	BMM=BMM+sl*sl.';
	end

 
 	logpdf = @(b1) +(n/2*log(det(D*sk(b1)*D*sigma0)))-1/2*sigma0*trace(sk((b1))*BMM)-vb*sum(1./b1)-(va+1)*sum(log((b1)))-sum(b1<=0)*10000000;
    rndd=slicesample(rndd.',1,'logpdf',logpdf,'burnin', 0);
    Psi=diag(((rndd)));

%%% Sign and scale adjustment for the PX model   
	signB=repmat(sign(diag(B(1:R,1:R))).',[V,1]);
	for indd=1:n
		signBB=diag(sign(diag(B(1:R,1:R))));
		aLambda(:,:,indd)=signBB*sqrt(Psi)*Lambda(:,:,indd)*sqrt(Psi)*signBB;
	end
	
	aaB=signB.*B*diag(1./sqrt(diag(Psi)));
	aGamma=Gamma*skron(signBB*sqrt(Psi),signBB*sqrt(Psi));
                                     
    for i=1:n  
        nL(:,:,i)=aaB*aLambda(:,:,i)*aaB';
    end
       
    Long=(L-nL);
	for indd=1:n
		LL(:,:,indd)=Long(:,:,indd)*Long(:,:,indd);
	end
	
	vv=trace1((LL));
    tLong=(L-nL);
	sigma=gamrnd(b1+n*V*(V+1)/4,1/(1/2*vv+b2));
	ssk=skron(Psi,Psi);
	for i=1:(R*(R+1)/2)
		sig_gam(i,i)=gamrnd(c1+p*R*(R+1)/4,1/(1/2*trace(Gamma'*Gamma)+c2));
	end
    
	if iter>burnin
        sL=sL+nL;
        sLambda=sLambda+aLambda;
        sGamma=sGamma+aGamma;
        sB=sB+aaB;
        devhat=devhat-2*(-1/2*sigma*vv-log(1/sigma*2*pi)*n*V*(V+1)/4);
        sigmahat=sigmahat+sigma;
    end
    
    toc;
end   

%%% Save Output into the structure 
    Output.nB=sB/(Niter-burnin);
    Output.nLambda=sLambda/(Niter-burnin);
    Output.nGamma=sGamma/(Niter-burnin);
    for j=1:n
		Output.nL=Output.nB*Output.nLambda(:,:,i)*(Output.nB)';
	end
	Output.sigmahat=sigmahat/(Niter-burnin);
    Output.err =  ERROR(Output.nLambda, nB, L, n);	% reconstruction errors
    toc;
end

%  Subroutine to calculate reconstruction errors                                                
function error = ERROR(Lambda, B, L, n)
	error = 0;
	A_const=cell(n);
	for i =1:n
		A_const{i} = B*Lambda(:,:,i)*B';
		TEMP = L(:,:,i) - A_const{i};
		error = error + norm(TEMP,'fro')/norm(L(:,:,i),'fro');
	end;
	error=error/n;
end

