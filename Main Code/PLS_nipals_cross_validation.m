function [output] = PLS_nipals_cross_validation(X, Y, translation)
% USAGE: [T,P,W,Wstar,U,b,C,Bpls,Bpls_star,Xhat,Yhat,R2X,R2Y]=PLS_nipals(X,Y,nfact)
% PLS regression NIPALS algorithm
% Y is actual set of outputs
% Y hat are the extimated values of the outputs

% COmpute the PLS regression coefficients
% X=T*P' Y=T*B*C'=X*Bpls  X and Y being Z-scores
%                          B=diag(b)
%    Y=X*Bpls_star with X being augmented with a col of ones
%                       and Y and X having their original units
% T'*T=I (NB normalization <> SAS)
% W'*W=I
%
% Test for PLS regression
% Herve Abdi November 2002/rev November 2004
%
%
% Version with T, W, and C being unit normalized
% U, P are not
% nfact=number of latent variables to keep
% default = rank(X)

    T = translation{1};
    P = translation{2};
    W = translation{3};
    Wstar = translation{4};
    U = translation{5};
    b = translation{6};
    C = translation{7};
    B_pls = translation{8};
    Bpls_star = translation{9};
    Xhat = translation{10};
    Yhat = translation{11};
    R2X = translation{12};
    R2Y = translation{13};

%     if nargin == 2
%         nfactor = rank(X);
%     end
%     X_ori=X;
%     Y_ori=Y;
%     if exist('nfactor')~=1;nfactor=rank(X);end
%     M_X=mean(X);
%     M_Y=mean(Y);
%     S_X=std(X);
%     S_Y=std(Y);
%     X=zscore(X);
%     Y=zscore(Y);
%     [nn,np]=size(X) ;
%     [n,nq]=size(Y)  ;
%     if nn~= n;
%         error(['Incompatible # of rows for X and Y']);
%     end
%     % Precision for convergence
%     epsilon=eps;
%     % # of components kepts
%     % Initialistion
%     % The Y set
%     U=zeros(n,nfactor);
%     C=zeros(nq,nfactor);
%     % The X set
%     T=zeros(n,nfactor);
%     P=zeros(np,nfactor);
%     W=zeros(np,nfactor);
%     b=zeros(1,nfactor);
%     R2_X=zeros(1,nfactor);
%     R2_Y=zeros(1,nfactor);
% 
%     Xres=X;
%     Yres=Y;
%     SS_X=sum(sum(X.^2));
%     SS_Y=sum(sum(Y.^2));
%     for l=1:nfactor
%         t=normaliz(Yres(:,1));
%         t0=normaliz(rand(n,1)*10);
%         u=t;
%         nstep=0;
%         maxstep=100;
%         while ( ( (t0-t)'*(t0-t) > epsilon/2) & (nstep < maxstep));
%             nstep=nstep+1;
%             %    disp(['Latent Variable #',int2str(l),'  Iteration #:',int2str(nstep)])
%             t0=t;
%             w=normaliz(Xres'*u);
%             t=normaliz(Xres*w);
%             % t=Xres*w;
%             c=normaliz(Yres'*t);
%             u=Yres*c;
%         end;
%         %% Take out display command for now
%         %   disp(['Latent Variable #',int2str(l),', convergence reached at step ',...
%         %     int2str(nstep)]);
% 
%         %X loadings
%         p=Xres'*t;
%         % b coef
%         b_l=((t'*t)^(-1))*(u'*t);
%         b_1=u'*t;
%         % Store in matrices
%         b(l)=b_l;
%         P(:,l)=p;
%         W(:,l)=w;
%         T(:,l)=t;
%         U(:,l)=u;
%         C(:,l)=c;
%         % deflation of X and Y
%         Xres=Xres-t*p';
%         Yres=Yres-(b(l)*(t*c'));
%         R2_X(l)=(t'*t)*(p'*p)./SS_X;
%         R2_Y(l)=(t'*t)*(b(l).^2)*(c'*c)./SS_Y;
% 
%     end

    
    Yhat=X*B_pls;
    X_rec=T*P';
    Y_rec=T*diag(b)*C';
%     Y_residual=Y-Y_rec;
    %% Bring back X and Y to their original units
    %
%     M_X=mean(X);
%     M_Y=mean(Y);
%     S_X=std(X);
%     S_Y=std(Y);
%     [n,~]=size(Y);
%     
%     Xori_rec=X_rec*diag(S_X)+(ones(n,1)*M_X);
%     Yori_rec=Y_rec*diag(S_Y)+(ones(n,1)*M_Y);
    %Unscaled_Y_hat=Yhat*diag(S_Y)+(ones(n,1)*M_Y);
    % The Wstart weights gives T=X*Wstar
    %
%     Wstar=W*inv(P'*W);
%     B_pls=Wstar*diag(b)*C';
    %% Bpls_star Y=[1|1|X]*Bpls_star
    %Bpls_star=[M_Y;[-M_X;eye(np,np)]*diag(S_X.^(-1))*B_pls*diag(S_Y)]
%     Bpls_star=[[-M_X;eye(np,np)]*diag(S_X.^(-1))*B_pls*diag(S_Y)];
%     Bpls_star(1,:)=Bpls_star(1,:)+M_Y;
%     
%     output = {T, P, W, Wstar, U, b, C, B_pls, Bpls_star, Xori_rec, Yori_rec, R2_X, R2_Y};
    output = {Yhat, X_rec, Y_rec};
end