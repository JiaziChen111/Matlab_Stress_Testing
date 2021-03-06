% /***************************************************
% **
% ** Procedures to analyze matrices:
% **
% ** x = dlyap(a,q)
% **
% */


function [ua,ta,ub,tb] = dlyap(a,c)
%/* Discrete Lyapunov Equation Solver
%** X = DLYAP(A,Q) solves the discrete Lyapunov Equation
%** AXA' - X+Q = 0
%*/

[ma,na] = size(a);
% ma = rows(a);
% na = cols(a);

[mc,nc] = size(c);
% mc = rows(c);
% nc = cols(c);

% /* Check dimension of matrices
% ** omitted
% */

a = inv( a+eye(ma) )*( a-eye(ma) );
c = 0.5*( eye(ma)-a )*c*( eye(ma)-a' );
b = a';

[mb,nb]=size(b);
% mb = rows(b);
% nb = cols(b);

% /* Perform Schur decomposition on A and convert to complex form
% */
[ua, ta] = schur(a);
[ub, tb] = schur(b);

% /* check all combinations of ta[i,i]+tb[j,j] for zero
% ** omitted
% */

% /* Transform C
% */
ucu = -ua'*c*ub;

% /* solve for first column of transformed solution
% */
y = zeros(ma,mb);
ema = eye(ma);
y(:,1) = inv(ta + ema*tb(1,1))*ucu(:,1);

% /* solve for remaining columns of transformed solution
% */
k = 2;
while k <= mb 
   km1 = 1:1:k-1;
   y(:,k) = inv( ta+ema*tb(k,k))*( ucu(:,k) - y(:,km1)*tb(km1,k));
   k = k+1;
end
   
% /* Find untransformed solution
% */
X = ua*y*ub';

% /* Ignore complex part
% */
X = real(X);

% /* Force X to be symmetric if C is symmetric
% */
if c == c';
   X = 0.5*(X+X');
end

