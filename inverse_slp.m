function [q0,nu0]=inverse_slp(a,b,m,n,L,LL,nue,lambda0,lambdaS,M1)%,A11, A21, B11,B21,A12 ,A22, B12, B22)
%% Barcilon's inverse SLP method using Magnus method for Sturm-Liouville problems of the form
% as described in the paper:
% Iterative solution of the inverse Sturm?Liouville problem
% Journal of Mathematical Physics 15, 429 (1974); https://doi.org/10.1063/1.1666664
% Victor Barcilon
% $$ (p(x)y')'+q(x)y=\lambda w(x) y $$
% p=w=1 and sets of BCs: y(0)=y(1)=0 and y'(0)=y(1)=0
%% Input:
%_a_ (double) left end point of the x-domain
%_b_ (double) right end point of the x-domain
%_lambdaS_ (double)  right end of eigenvalue interval
%_lambda0_ (double) left end of eigenvalue interval
%_m_ (integer) # divisions for the eigenvalue interval
%_n_ (integer) # divisions for the [a,b]
%_L_ (integer) #  multisection iteration steps
%_LL_ (integer) #  length of exact eigenvalue sequence
%_M1_ (integer) #  inverses iteration steps
%_nue_ (array) #  exact eigenvalue sequence
%% Output:
%_q0_ (function) reconstructed potential
%_nu0_ (array) reconstructed eigenvalue sequence
% Author: Upeksha Perera (April 2019)
% Supplementary Material for the article titled
% "Solutions of Direct and Inverse Even-order Sturm-Liouville problems using Magnus expansion"
% by Upeksha Perera and Christine Böckmann
% Correspondence: upeksha@kln.ac.lk; Department of Mathematics, University of Kelaniya, 11600 Kelaniya, Sri Lanka;
% Current address: Institut für Mathematik, Universität Potsdam, 14476 Potsdam, Germany; bodhiyabadug@uni-potsdam.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% USAGE:
% % Example: Reconstructing the potential for  the  problem
% % -y''+cos(t)*y=lambda y,
% % given the corresponding eigenvalue sequences for y(0)=y(1)=0 and y'(0)=y(1)=0
% %
% a=0; b=pi; % end points
% lambda0=0; lambdaS=20;   % eigenvalue search interval
% m=10;n=10; L=2; % parameters for magnus method
% A11=1; A21= 0; B11=1 ; B21=0 ; % first set of BCs
% A12=0  ; A22= 1; B12= 1 ; B22=0; % second set of BCs
% q0e=@(t) cos(t); % exact potential
% p=@(t) -ones(size(t)); w0=@(t) ones(size(t)); % these two coefficients are fixed
% M1=2; % number of inverse algorithm steps

% lambdae=slpm(a,b,lambda0,lambdaS,p,q0e,w0,m,n,L,A11,A21,B11,B21);
% mue=slpm(a,b,lambda0,lambdaS,p,q0e,w0,m,n,L,A12,A22,B12,B22);
% LL=min([length(lambdae),length(mue)]);
% for k=1:LL
%     nue(2*k-1,:)=4*mue(k);
%     nue(2*k,:)=4*lambdae(k);
% end
% % nue: exact eigenvalues
% 
% [q0,nu0]=inverse_slp(a,b,m,n,L,LL,nue,lambda0,lambdaS,M1);%,A11, A21, B11,B21,A12 ,A22, B12, B22);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A11=1; A21= 0; B11=1 ; B21=0 ; % first set of BCs
A12=0  ; A22= 1; B12= 1 ; B22=0; % second set of BCs

x=linspace(a,b,n); % partitioning x-domain
sum1=zeros(M1,length(x)); % store the error terms

%% constructing w and omega for q=0
for k=1:2*lambdaS
    cc1=sqrt(2/(3*b));
    w(k,:)=(cos(k.*x*pi/b)-cos(k*pi))*cc1;
    dw(k,:)=-(k*pi/b)*sin(k.*x*pi/b)*cc1;
    
    cc2=sqrt(2/b);
    omega(k,:)=sin(k.*x*pi/b)*cc2;
    domega(k,:)=(k*pi/b)*cos(k.*x*pi/b)*cc2;
    
    int1(k)=-sqrt(2*b/3)*cos(k*pi);
    int2(k)=k*pi/(b*sqrt(3));
    c(k)=int1(k)./int2(k);
    dw1(k,:)=c(k)*domega(k,:);
end

p=@(t) -ones(size(t)); w0=@(t) ones(size(t)); % these two coefficients are fixed
q0=@(t) zeros(size(t));% initial guess

lambda=slpm(a,b,lambda0,lambdaS,p,q0,w0,m,n,L,A11,A21,B11,B21);
mu=slpm(a,b,lambda0,lambdaS,p,q0,w0,m,n,L,A12,A22,B12,B22);

for k=1:LL
    nu0(2*k-1,:)=4*mu(k);
    nu0(2*k,:)=4*lambda(k);
end
% nu0: approximate eigenvalues
%%
for M=1:M1 % iteration loop
    
    for kk=1:LL
        sum1(M,:)=sum1(M,:)+dw1(kk,:).*(nue(kk)-nu0(kk));
    end
    sum1(M,:)=0.25*sum1(M,:);
    
    Q0=q0(x)+sum1(M,:); % updating q
    pp=griddedInterpolant(x',Q0','pchip'); % getting the functional form of q
    q0=@(t) pp(t);
    
    % calculating new eigenvalues using updated q
    lambda=[]; mu=[]; nu0=[];
    lambda=slpm(a,b,lambda0,lambdaS,p,q0,w0,m,n,L,A11,A21,B11,B21);
    mu=slpm(a,b,lambda0,lambdaS,p,q0,w0,m,n,L,A12,A22,B12,B22);
    for k=1:LL
        nu0(2*k-1,:)=4*mu(k);
        nu0(2*k,:)=4*lambda(k);
    end
end
