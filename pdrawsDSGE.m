function ymat = pdrawsDSGE(para,state,shocks,steps)

[TT,QQ,RR,HH,DD,ZZ,VV,RC]=sysmat(para);

s = state';

ymat = [];

shocks = shocks * chol(QQ);

for i = 1:steps
    
    sn  = TT * s + RR * shocks(i,:)';
    
    ymat(:,i)   = DD + ZZ * sn;
  
    s = sn;
    %s(16) = 0.001;
    
end