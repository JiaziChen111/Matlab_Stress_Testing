function seq=seqa(a,b,c);
% PURPOSE: produce a sequence of values
% -----------------------------------------------------
% USAGE: y = seqa(a,b,c)
%  where    a = initial value in sequence 
%           b = increment
%           c = number of values in the sequence  
% -----------------------------------------------------
% RETURNS: a sequence, (a:b:(a+b*(c-1)))' in MATLAB notation
% ----------------------------------------------------- 
% NOTE: a Gauss compatability function
% -----------------------------------------------------

% written by:
% James P. LeSage, Dept of Economics
% Texas State University-San Marcos
% 601 University Drive
% San Marcos, TX 78666
% jlesage@spatial-econometrics.com
       
% seqa Gauss eqivalent of seqa(a,b,c)
seq=(a:b:(a+b*(c-1)))';
return;
