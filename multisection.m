function y=multisection(f,a,b,tol,m)
%% Multisection method for finding a root in an interval
% Input:
%_f_ (function handle)
%_a_ (double) left end point of the interval
%_b_ (double) right end point of the interval
%_tol_ (double) tolerance
%_m_ (integer) # maximum iterations
% Output:
%_y_ (double) approximate root of $f$ in the interval [a, b] with a tolerance of 'tol'
% and maximum number of iterations m
% Author: Upeksha Perera (April 2019)
% Supplementary Material for the article titled
% "Solutions of Direct and Inverse Even-order Sturm-Liouville problems using Magnus expansion"
% by Upeksha Perera and Christine Böckmann
% Correspondence: upeksha@kln.ac.lk; Department of Mathematics, University of Kelaniya, 11600 Kelaniya, Sri Lanka;
% Current address: Institut für Mathematik, Universität Potsdam, 14476 Potsdam, Germany; bodhiyabadug@uni-potsdam.de
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % USAGE:
% % Example: Finding the  root  in the  interval [1,2]  for  the  function $f=x^3-x-2$
% % with a tolerance 1e-7 and maximum  of 10 iterations
% % a=1;b=2;
% % f=@(x) x.^3-x-2;
% % tol=1e-7;
% % m=10;
% % y=multisection(f,a,b,tol,m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=(b-a)/m;   % step -size
x=a:h:b; % sequence  of  numbers: x1=a , x2=a+h , ... , xm=b

values=zeros(1,m);  % store function values at each subdivision point

err = 1;
while err > tol
    for k=1:m+1
        values(k)=f(x(k));  % calculate function value at each point k
    end
    
    scinter = find(diff(sign(values))); % find the sign changing position of f
    
    % reset the right and left end points to bracket the sign changing position
    b=x(scinter(1)+1);
    a=x(scinter(1));
    
    % redefine the step size
    h=(b-a)/m;
    
    % refine the root interval
    x=a:h:b;
    
    err=abs(f(x));
end % end of while loop (multisection)
y=(a+b)/2; % approximate root