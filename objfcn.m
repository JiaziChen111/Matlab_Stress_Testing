function fx = objfcn(para,v)

% Inputted parameter values are transformations of model parameters for
% puroposes of minimizing the obj function. Here, we transform them back
% to model values for inputting into evalmod
modelpara = trans(para,v.trspec);

% If any parameters are fixed, they are set here:
modelpara = modelpara.*v.pmaskinv + v.pfix.*v.pmask;
% 
[lnpy,retcode,obserror,obsvar] = evalmod(modelpara,v.data);

% 
 if isnan(lnpy);
    lnpy = -1E8;
 end;
% 
% /** Evaluate the Prior distribution **/
 lnprio = priodens(modelpara, v.pmean, v.pstdd, v.pshape);
% 
% /** minimize the negative of log posterior **/
fx = real(-lnpy-lnprio);

