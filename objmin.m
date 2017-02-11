function[p,g,csretcode] = objmin(p,v)

% Use Sims' csminwel to maximize the log likelihood function 'objfcn'.
% Input parameters are:
%
%   p0vals:     Starting values for maximizer
%   p:          A structure that holds starting values of the parameters
%               and estimated values of the parameters
%   v:          A structure of parameters/values that get passed to 
%               csminwel

npara  = length(p);
p0vals = zeros(npara,1);
cc     = 0.1;

for i = 1:npara
    p0vals(i) = p(i).sval;
end

% Use invtrans to change the parameter support vector so minimization
% routine is less likely to hang. x0 is the resulting start val vector

x0 = invtrans(p0vals,v.trspec) + cc* ( randn(npara,1) );

% Initial Hessian
H0 = eye(npara)*1E-4;

% Optimization environment parameters
nit  = 2000;        
crit = 1E-8;        

% call csminwel to maximize the log likelihood
[fh,xh,g,H,itct,fcount,csretcode] = csminwel('objfcn',x0,H0,[],crit,nit,v);
%  from max back to model parameters using trans
paraest = trans(xh,v.trspec);
paraest = real(paraest);
paraest = paraest.*v.pmaskinv + v.pfix.*v.pmask;

for i = 1:npara
    p(i).estval = paraest(i);
end
