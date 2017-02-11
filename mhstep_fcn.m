%% MH sampling function definition

function [lnpost,lnpy] = mhstep_fcn(para,v)

bounds = v.trspec(:,[2,3]);

ind1 = para > bounds(:,1);
ind2 = para < bounds(:,2);

ind1 = delif(ind1,v.pmask);
ind2 = delif(ind2,v.pmask);

if all(ind1) && all(ind2)
    
    [lnpy,retcode,obserror,obsvar] = evalmod(para,v.data);
    lnpy = real(lnpy);
    lnprio = priodens(para, v.pmean, v.pstdd, v.pshape);
    lnpost = lnpy + real(lnprio);
    
else
    
    lnpost = -1e6;
    lnpy = -1e20;
     
end