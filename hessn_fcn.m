%% Objective function for hessian computation

function fx = hessn_fcn(para,v)

%modelpara = trans(para,v.trspec);
modelpara = para.*v.pmaskinv + v.pfix.*v.pmask;

[lnpy,retcode,obserror,obsvar] = evalmod(modelpara,v.data);

lnprio = priodens(modelpara,v.pmean,v.pstdd,v.pshape);

fx = lnpy + lnprio;

