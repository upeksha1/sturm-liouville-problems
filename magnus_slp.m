function ev=magnus_slp(a,b,lambda0,lambdaS,p,q,w,m,n,LM,A1,A2,B1,B2)

%% Magnus method for Sturm-Liouville problems of the form
% $$ (p(x)y')'+q(x)y=\lambda w(x) y $$
% Input:
%_a_ (double) left end point of the x-domain
%_b_ (double) right end point of the x-domain
%_lambdaS_ (double)  right end of eigenvalue interval
%_lambda0_ (double) left end of eigenvalue interval
%_q_ (function) coefficient of y
%_p_ (function) coefficient of y''
%_w_ (function) coefficient of lambda y term
%_m_ (integer) # divisions for the eigenvalue interval
%_n_ (integer) # divisions for the [a,b]
%_LM_ (integer) #  multisection iteration steps
%_A1_ (0 or 1) coefficient of boundary condition y(a)
%_A2_ (0 or 1) coefficient of boundary condition y'(a)
%_B1_ (0 or 1) coefficient of boundary condition y(b)
%_B2_ (0 or 1) coefficient of boundary condition y'(b)
%% note that A1,A2,B1,B2  selected such that
% A1*A2'=A2*A1', B1*B2'=B2*B1' and (A1,A2) and (B1,B2) have rank 1
% Output:
%_ev_ (double) approximate eigenvalue in the interval [lambda0, lambdaS]
% Author: Upeksha Perera (April 2019)
% Supplementary Material for the article titled
% "Solutions of Direct and Inverse Even-order Sturm-Liouville problems using Magnus expansion"
% by Upeksha Perera and Christine Böckmann
% Correspondence: upeksha@kln.ac.lk; Department of Mathematics, University of Kelaniya, 11600 Kelaniya, Sri Lanka;
% Current address: Institut für Mathematik, Universität Potsdam, 14476 Potsdam, Germany; bodhiyabadug@uni-potsdam.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % USAGE:
% % Example 1: Finding the  eigenvalue  in the  interval [0, 10]  for  the  problem
% % y''+lambda y=0, y(0)=y(1)=0
%
% % a=0; b=1;
% % lambda0=0; lambdaS=10;
% % q=@(x) 0;
% % p=@(x) 1;
% % w=@(x) 1;
% % m=100; n=100;
% % LM=5;
% % A1=1; A2=0;
% % B1=1; B2=0;
% % ev=magnus_slp(a,b,lambda0,lambdaS,p,q,w,m,n,LM,A1,A2,B1,B2);

% % Example 2: Paine problem 2 in Pryce [1993] from Paine et al. [1981].
% % y''+1/(x+0.1)^2 y=lambda y, y(0)=y(pi)=0
% a=0; b=pi;
% lambda0=0; lambdaS=2;
% q=@(x) 1./(x+0.1).^2;
% p=@(x) -1;
% w=@(x) 1;
% m=10; n=10; L=5;
% A1=1; A2=0;
% B1=1; B2=0;
% ev=magnus_slp(a,b,lambda0,lambdaS,p,q,w,m,n,L,A1,A2,B1,B2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h1=(lambdaS-lambda0)/m;  % step-size of eigenvalue interval
lambda=lambda0:h1:lambdaS; % sequence of numbers: lambda0, lambda1, ... , lambdaS
h=(b-a)/n; % step-size for the domain
x=a:h:b; % sequence of numbers a=x0, x1, ..., xn=b

Y=eye(2);  % initial value of Y(0)

A=[A1,A2; 0 0];
B=[0 0;B1,B2];

G=@(x,lambda)[0  1./p(x) ;  lambda.*w(x)-q(x) 0]; % matrix G

c1=1/2-sqrt(15)/10; % Gaussian point c1
c3=1/2+sqrt(15)/10; % Gaussian point c3

G1=@(x,lambda) G(x+c1.*h,lambda); % G(c1*h)
G2=@(x,lambda) G(x+h/2,lambda); % G(c2*h), c2=1/2
G3=@(x,lambda) G(x+c3.*h,lambda); % G(c3*h)

Q1=@(x,lambda) h.*G2(x,lambda); % Q1=h*G(c1*h)
Q2=@(x,lambda) (sqrt(15)*h/3).*(G3(x,lambda)-G1(x,lambda));
Q3=@(x,lambda) (10*h/3).*(G3(x,lambda)-2*G2(x,lambda)+G1(x,lambda));

R1=@(x,lambda)Q1(x,lambda)*Q2(x,lambda)-Q2(x,lambda)*Q1(x,lambda);
% R1=[Q1,Q2]=Q1*Q2-Q2*Q1

R2=@(x,lambda) Q1(x,lambda)*(2.*Q3(x,lambda)+R1(x,lambda)) ...
    -(2.*Q3(x,lambda)+R1(x,lambda))*Q1(x,lambda);  % R2=[Q1,2Q3+R1]

R3=@(x,lambda)(-20*Q1(x,lambda)-Q3(x,lambda)+R1(x,lambda))*(Q2(x,lambda) ...
    -R2(x,lambda)./60)-(Q2(x,lambda)-R2(x,lambda)./60)*(-20*Q1(x,lambda) ...
    -Q3(x,lambda)+R1(x,lambda));

sigma=@(x,lambda) Q1(x,lambda)+Q3(x,lambda)./12+R3(x,lambda)./240;

%%%%%%%%%%%% the multisection method starts %%%%%%%%%
L=0;
while L<LM % # multisection steps
    for k=1:m+1
        Y=eye(2);
        for j=1:n
            Y=expm(sigma(x(j),lambda(k)))*Y; % Calculate Y(b) iteratively
        end
        F(k)=det(A+B*Y);  % calculate characteristic function at each point k
    end
    
    scinter = find(diff(sign(F))); % find the sign changing position of F
    
    % reset the right and left end points to bracket the sign changing position
    lambdaS=lambda(scinter(1)+1);
    lambda0=lambda(scinter(1));
    
    % redefine the step size
    h1=(lambdaS-lambda0)/m;
    
    % refine the possible lambda values
    lambda=lambda0:h1:lambdaS;
    
    L=L+1; % increase the multisection counter
end % end of while loop (multisection)
%%%%%%%%%%%%%%%%% end of multisection method %%%%%%%%%%

ev=(lambda0+lambdaS)/2; % approximate eigenvalue